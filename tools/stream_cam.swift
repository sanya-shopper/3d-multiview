// multiview live camera streamer: webcam -> grayscale -> TCP to livehub.
// Runs on each camera MacBook; macOS system frameworks only.
//
//   swiftc -O tools/stream_cam.swift -o stream_cam
//   ./stream_cam <hub-host> <port> <camid> [fps]
//
// First run prompts for camera permission (grant it in System Settings
// if the prompt is missed). Status changes are also SPOKEN (macOS say)
// because mid-session this screen faces the volume, not the operator;
// MV_MUTE=1 silences. Frames go out as:
//   "MVFR" | u32 camid | u32 w | u32 h | u64 seq | f64 t_mono | w*h gray
// little-endian, t_mono = this machine's monotonic uptime seconds (each
// camera is its own clock domain; the display counter is the shared
// clock, per the paper's temporal model).
import AVFoundation
import AppKit
import CoreGraphics
import Foundation

// Every argument optional with deployment defaults:
//   secondlaptop [hub-ip] [port] [camid] [fps]   (172.20.10.2 9900 2 5)
let a = CommandLine.arguments
let host = a.count > 1 ? a[1] : "172.20.10.2"
let port = a.count > 2 ? (UInt16(a[2]) ?? 9900) : 9900
let camid = a.count > 3 ? (UInt32(a[3]) ?? 2) : 2
let fps = a.count > 4 ? (Double(a[4]) ?? 5.0) : 5.0

// Spoken feedback: fire-and-forget /usr/bin/say child; it survives an
// exit() right after, so failure paths can speak then die.
func speak(_ text: String) {
    if ProcessInfo.processInfo.environment["MV_MUTE"] != nil { return }
    let p = Process()
    p.executableURL = URL(fileURLWithPath: "/usr/bin/say")
    p.arguments = [text]
    try? p.run()
}

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
// One lock serializes every message onto the socket: frames go out on
// the capture queue while clock-sync replies (MVTS, below) go out on
// the main-thread timer, and interleaved partial send()s would tear a
// frame mid-stream (the hub resyncs on magic, but the frame is lost).
let sendLock = NSLock()
func sendAll(_ s: Int32, _ data: [UInt8]) -> Bool {
    sendLock.lock()
    defer { sendLock.unlock() }
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
// Push a diagnostics line to the hub over the same socket (MVLG
// message): the hub writes it to rec/camlog_<id>.txt so all logs end up
// on the hub machine with no SSH or file transfer.
func sendLog(_ s: Int32, _ camid: UInt32, _ text: String) {
    let body = Array(text.utf8)
    var msg: [UInt8] = [0x4D, 0x56, 0x4C, 0x47] // "MVLG"
    func le32(_ v: UInt32) {
        msg.append(UInt8(v & 255)); msg.append(UInt8((v >> 8) & 255))
        msg.append(UInt8((v >> 16) & 255)); msg.append(UInt8((v >> 24) & 255))
    }
    le32(camid); le32(UInt32(body.count))
    msg.append(contentsOf: body)
    _ = sendAll(s, msg)
}

// retry for ~15 s so the second laptop can be started before OR after
// the hub without failing
for _ in 0..<150 {
    sock = connectHub(host, port)
    if sock >= 0 { break }
    usleep(100_000)
}
guard sock >= 0 else {
    print("cannot reach hub at \(host):\(port) -- check the IP, same Wi-Fi, hub running")
    speak("camera \(camid) cannot reach the hub")
    exit(1)
}
print("connected to \(host):\(port) as camera \(camid), \(fps) fps")

// read hub messages from the socket: MVAK acks (rx | decoded, so the
// sender knows the hub is really receiving, not just that TCP accepted
// the bytes) and MVPB clock-sync probes ("MVPB" | f64 t_hub_send),
// answered immediately with "MVTS" | u32 camid | f64 t_hub_echo | f64
// t_cam.  Both hub->camera messages are 12 bytes but can arrive SPLIT
// at any byte boundary, so a persistent carry buffer (g.rxPending)
// holds the partial tail between timer ticks; unrecognized bytes are
// skipped one at a time (resync on the next magic).
func readAcks(_ g: Grabber) {
    var buf = [UInt8](repeating: 0, count: 4096)
    let n = buf.withUnsafeMutableBytes { recv(sock, $0.baseAddress, 4096, MSG_DONTWAIT) }
    if n > 0 { g.rxPending.append(contentsOf: buf[0..<n]) }
    let bytes = g.rxPending
    var i = 0
    while i + 12 <= bytes.count {
        if bytes[i] == 0x4D && bytes[i+1] == 0x56 && bytes[i+2] == 0x41 && bytes[i+3] == 0x4B {
            func u32(_ o: Int) -> UInt32 {
                UInt32(bytes[i+o]) | (UInt32(bytes[i+o+1])<<8)
                | (UInt32(bytes[i+o+2])<<16) | (UInt32(bytes[i+o+3])<<24) }
            g.hubRx = u32(4); g.hubDecoded = u32(8)
            g.lastAck = ProcessInfo.processInfo.systemUptime
            i += 12
        } else if bytes[i] == 0x4D && bytes[i+1] == 0x56 && bytes[i+2] == 0x50 && bytes[i+3] == 0x42 {
            // MVPB: echo t_hub_send BIT-EXACTLY (raw byte copy, never
            // through a Double) and stamp t_cam at reply time; the
            // hub midpoints the probe, so reply latency only widens
            // RTT a little while a corrupted echo would poison the fit
            var reply: [UInt8] = [0x4D, 0x56, 0x54, 0x53] // "MVTS"
            var v = camid
            for _ in 0..<4 { reply.append(UInt8(v & 255)); v >>= 8 }
            reply.append(contentsOf: bytes[(i+4)..<(i+12)])
            var t = ProcessInfo.processInfo.systemUptime.bitPattern
            for _ in 0..<8 { reply.append(UInt8(t & 255)); t >>= 8 }
            _ = sendAll(sock, reply)
            i += 12
        } else { i += 1 }
    }
    g.rxPending.removeFirst(i) // keep the (< 12 byte) partial tail
}

// ---- capture ----
final class Grabber: NSObject, AVCaptureVideoDataOutputSampleBufferDelegate {
    var seq: UInt64 = 0
    var lastSent = 0.0
    var status = "starting"
    var dims = "?"
    var hubRx: UInt32 = 0        // frames the hub says it received
    var hubDecoded: UInt32 = 0
    var lastAck = 0.0           // uptime of last ACK from the hub
    var rxPending: [UInt8] = [] // partial hub message carry (readAcks)
    var sendOK = true
    let minGap: Double
    let camid: UInt32
    let sock: Int32
    init(minGap: Double, camid: UInt32, sock: Int32) {
        self.minGap = minGap
        self.camid = camid
        self.sock = sock
    }
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
        le32(self.camid); le32(UInt32(w)); le32(UInt32(h))
        var s = seq
        for _ in 0..<8 { msg.append(UInt8(s & 255)); s >>= 8 }
        var t = now.bitPattern
        for _ in 0..<8 { msg.append(UInt8(t & 255)); t >>= 8 }
        msg.append(contentsOf: gray)
        if !sendAll(self.sock, msg) {
            sendOK = false
            status = "HUB CONNECTION LOST"
            print("hub connection lost")
            speak("camera \(self.camid) lost the hub")
            exit(1)
        }
        sendOK = true
        seq += 1
        dims = "\(w)x\(h)"
        status = "streaming to hub"
        if seq == 1 {
            sendLog(sock, camid, "camera \(camid) connected, \(w)x\(h)")
        }
        if seq % 50 == 0 {
            print("sent \(seq) frames (\(w)x\(h))")
            sendLog(sock, camid, "sent \(seq) frames, \(w)x\(h)")
        }
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
let grabber = Grabber(minGap: 1.0 / fps, camid: camid, sock: sock)
out.setSampleBufferDelegate(grabber,
                            queue: DispatchQueue(label: "grab"))
session.addOutput(out)
session.startRunning()
print("streaming from \(dev.localizedName); close the window or Ctrl-C to stop")

// ---- live preview window so the operator can see what the camera
//      sees and aim it correctly ----
let app = NSApplication.shared
app.setActivationPolicy(.regular)
let win = NSWindow(contentRect: NSRect(x: 80, y: 80, width: 720, height: 460),
                   styleMask: [.titled, .closable, .resizable, .miniaturizable],
                   backing: .buffered, defer: false)
win.title = "multiview camera \(camid)  ->  \(host):\(port)"
let preview = AVCaptureVideoPreviewLayer(session: session)
preview.videoGravity = .resizeAspect
let root = NSView(frame: win.contentRect(forFrameRect: win.frame))
root.wantsLayer = true
root.layer = preview
root.autoresizingMask = [.width, .height]
// status banner along the bottom
let banner = NSTextField(labelWithString: "connecting…")
banner.frame = NSRect(x: 0, y: 0, width: 720, height: 28)
banner.alignment = .center
banner.font = NSFont.monospacedSystemFont(ofSize: 13, weight: .bold)
banner.textColor = .white
banner.backgroundColor = NSColor.black.withAlphaComponent(0.65)
banner.drawsBackground = true
banner.autoresizingMask = [.width, .maxYMargin]
root.addSubview(banner)
win.contentView = root
win.makeKeyAndOrderFront(nil)
app.activate(ignoringOtherApps: true)

// banner: GREEN = hub is acknowledging our frames; AMBER = waiting for
// the camera / first frames; RED = sending but the hub isn't confirming
// State TRANSITIONS are spoken (once each), so the operator can aim and
// troubleshoot without seeing this screen.
var saidLive = false, saidDecoding = false, saidPermission = false
var inRed = false
let tStart = ProcessInfo.processInfo.systemUptime
Timer.scheduledTimer(withTimeInterval: 0.5, repeats: true) { _ in
    readAcks(grabber)
    let nowU = ProcessInfo.processInfo.systemUptime
    let ackAge = grabber.lastAck > 0 ? nowU - grabber.lastAck : 999
    if grabber.seq == 0 {
        banner.textColor = .systemYellow
        banner.stringValue = "waiting for camera… (grant camera permission if prompted)"
        if !saidPermission && nowU - tStart > 5 {
            saidPermission = true
            speak("camera \(camid) is waiting for camera permission")
        }
    } else if !grabber.sendOK || ackAge > 4 {
        banner.textColor = .systemRed
        banner.stringValue =
            "⚠ HUB NOT CONFIRMING  camera \(camid)  ·  sent \(grabber.seq)"
            + "  ·  check the hub is running and on this Wi-Fi"
        if !inRed {
            inRed = true
            speak("camera \(camid), hub not confirming")
        }
    } else {
        banner.textColor = .systemGreen
        banner.stringValue =
            "● LIVE  camera \(camid)  \(grabber.dims)  ·  sent \(grabber.seq)"
            + "  ·  hub got \(grabber.hubRx), decoded \(grabber.hubDecoded)"
            + "  ·  aim so the pattern fills the view"
        if inRed {
            inRed = false
            saidLive = true
            speak("camera \(camid) confirmed again")
        } else if !saidLive {
            saidLive = true
            speak("camera \(camid) live, hub confirming")
        }
        // the aiming cue: the hub decoding frames means the pattern is
        // actually in this camera's view
        if !saidDecoding && grabber.hubDecoded > 0 {
            saidDecoding = true
            speak("hub is decoding camera \(camid)")
        }
    }
}
app.run()
