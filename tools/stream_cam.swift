// multiview live camera streamer: webcam -> grayscale -> TCP to livehub.
// Runs on each camera MacBook; macOS system frameworks only.
//
//   swiftc -O tools/stream_cam.swift -o stream_cam
//   ./stream_cam <hub-host> <port> <camid> [fps]
//
// First run prompts for camera permission (grant it in System Settings
// if the prompt is missed). Frames go out as:
//   "MVFR" | u32 camid | u32 w | u32 h | u64 seq | f64 t_mono | w*h gray
// little-endian, t_mono = this machine's monotonic uptime seconds (each
// camera is its own clock domain; the display counter is the shared
// clock, per the paper's temporal model).
import AVFoundation
import CoreGraphics
import Foundation

let a = CommandLine.arguments
guard a.count == 4 || a.count == 5, let port = UInt16(a[2]),
      let camid = UInt32(a[3]) else {
    print("usage: stream_cam <hub-host> <port> <camid> [fps]")
    exit(1)
}
let fps = a.count == 5 ? (Double(a[4]) ?? 5.0) : 5.0

// ---- TCP client (BSD sockets) ----
var sock: Int32 = -1
func connectHub(_ host: String, _ port: UInt16) -> Int32 {
    var hints = addrinfo(ai_flags: 0, ai_family: AF_INET,
                         ai_socktype: SOCK_STREAM, ai_protocol: 0,
                         ai_addrlen: 0, ai_canonname: nil, ai_addr: nil,
                         ai_next: nil)
    var res: UnsafeMutablePointer<addrinfo>?
    guard getaddrinfo(host, String(port), &hints, &res) == 0,
          let r = res else { return -1 }
    defer { freeaddrinfo(res) }
    let s = socket(r.pointee.ai_family, r.pointee.ai_socktype,
                   r.pointee.ai_protocol)
    if s < 0 { return -1 }
    if connect(s, r.pointee.ai_addr, r.pointee.ai_addrlen) != 0 {
        close(s)
        return -1
    }
    var one: Int32 = 1
    setsockopt(s, IPPROTO_TCP, TCP_NODELAY, &one,
               socklen_t(MemoryLayout<Int32>.size))
    return s
}
func sendAll(_ s: Int32, _ data: [UInt8]) -> Bool {
    var off = 0
    while off < data.count {
        let n = data.withUnsafeBytes { p -> Int in
            send(s, p.baseAddress! + off, data.count - off, 0)
        }
        if n <= 0 { return false }
        off += n
    }
    return true
}

sock = connectHub(a[1], port)
guard sock >= 0 else {
    print("cannot connect to \(a[1]):\(port)")
    exit(1)
}
print("connected to \(a[1]):\(port) as camera \(camid), \(fps) fps")

// ---- capture ----
final class Grabber: NSObject, AVCaptureVideoDataOutputSampleBufferDelegate {
    var seq: UInt64 = 0
    var lastSent = 0.0
    let minGap: Double
    init(minGap: Double) { self.minGap = minGap }
    func captureOutput(_ output: AVCaptureOutput,
                       didOutput sampleBuffer: CMSampleBuffer,
                       from connection: AVCaptureConnection) {
        let now = ProcessInfo.processInfo.systemUptime
        if now - lastSent < minGap { return }
        lastSent = now
        guard let pb = CMSampleBufferGetImageBuffer(sampleBuffer) else {
            return
        }
        CVPixelBufferLockBaseAddress(pb, .readOnly)
        defer { CVPixelBufferUnlockBaseAddress(pb, .readOnly) }
        let w = CVPixelBufferGetWidth(pb), h = CVPixelBufferGetHeight(pb)
        let stride = CVPixelBufferGetBytesPerRow(pb)
        guard let base = CVPixelBufferGetBaseAddress(pb) else { return }
        // BGRA -> gray (Rec.601 integer approximation)
        var gray = [UInt8](repeating: 0, count: w * h)
        let src = base.assumingMemoryBound(to: UInt8.self)
        for y in 0..<h {
            let row = src + y * stride
            for x in 0..<w {
                let b = Int(row[4 * x]), g = Int(row[4 * x + 1]),
                    r = Int(row[4 * x + 2])
                gray[y * w + x] =
                    UInt8((77 * r + 150 * g + 29 * b) >> 8)
            }
        }
        var msg = [UInt8]()
        msg.reserveCapacity(32 + w * h)
        msg.append(contentsOf: [0x4D, 0x56, 0x46, 0x52]) // "MVFR"
        func le32(_ v: UInt32) {
            msg.append(UInt8(v & 255)); msg.append(UInt8((v >> 8) & 255))
            msg.append(UInt8((v >> 16) & 255))
            msg.append(UInt8((v >> 24) & 255))
        }
        le32(camid); le32(UInt32(w)); le32(UInt32(h))
        var s = seq
        for _ in 0..<8 { msg.append(UInt8(s & 255)); s >>= 8 }
        var t = now.bitPattern
        for _ in 0..<8 { msg.append(UInt8(t & 255)); t >>= 8 }
        msg.append(contentsOf: gray)
        if !sendAll(sock, msg) {
            print("hub connection lost")
            exit(1)
        }
        seq += 1
        if seq % 50 == 0 { print("sent \(seq) frames (\(w)x\(h))") }
    }
}

let session = AVCaptureSession()
session.sessionPreset = .hd1280x720
guard let dev = AVCaptureDevice.default(for: .video),
      let input = try? AVCaptureDeviceInput(device: dev) else {
    print("no camera (check permission)")
    exit(1)
}
session.addInput(input)
let out = AVCaptureVideoDataOutput()
out.videoSettings =
    [kCVPixelBufferPixelFormatTypeKey as String: kCVPixelFormatType_32BGRA]
let grabber = Grabber(minGap: 1.0 / fps)
out.setSampleBufferDelegate(grabber,
                            queue: DispatchQueue(label: "grab"))
session.addOutput(out)
session.startRunning()
print("streaming from \(dev.localizedName); Ctrl-C to stop")
RunLoop.main.run()
