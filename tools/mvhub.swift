// multiview HUB, one binary: launches the processing hub (livehub) as a
// child, captures its output into the collection AND onto the screen,
// and runs THIS Mac as camera 1 with a live preview window. The hub
// operator runs one thing.
//
//   swiftc -O tools/mvhub.swift -o mvhub     (needs ./livehub beside it)
//   ./mvhub [port] [pitch_mm] [record_dir] [fps]
// defaults: 9900 0.1133 rec 5
import AVFoundation
import AppKit
import Foundation

let a = CommandLine.arguments
let port = a.count > 1 ? (UInt16(a[1]) ?? 9900) : 9900
let pitch = a.count > 2 ? a[2] : "0.1133"
let recdir = a.count > 3 ? a[3] : "rec"
let fps = a.count > 4 ? (Double(a[4]) ?? 5.0) : 5.0
let camid: UInt32 = 1
let kitDir = URL(fileURLWithPath: a[0]).deletingLastPathComponent()

// ---- launch the hub child, tee its output to terminal + rec/hub.log ----
try? FileManager.default.createDirectory(atPath: recdir,
                                         withIntermediateDirectories: true)
let logPath = "\(recdir)/hub.log"
FileManager.default.createFile(atPath: logPath, contents: nil)
let logfh = try? FileHandle(forWritingTo: URL(fileURLWithPath: logPath))
let hub = Process()
hub.executableURL = kitDir.appendingPathComponent("hubengine")
hub.arguments = [String(port), pitch, recdir]
let pipe = Pipe()
hub.standardOutput = pipe
hub.standardError = pipe
pipe.fileHandleForReading.readabilityHandler = { h in
    let d = h.availableData
    if d.isEmpty { return }
    FileHandle.standardOutput.write(d)   // to the terminal
    logfh?.write(d)                      // into the collection
}
do { try hub.run() } catch {
    FileHandle.standardError.write(Data("cannot start ./hubengine (must be beside hublaptop)\n".utf8))
    exit(1)
}

// Spoken feedback for THIS process's own failure paths; the hubengine
// child speaks the hub-level events (ready/connected/calibrated/paired)
// itself on this same Mac, so only camera-1 local issues live here.
// MV_MUTE=1 silences (the child inherits the environment).
func speak(_ text: String) {
    if ProcessInfo.processInfo.environment["MV_MUTE"] != nil { return }
    let p = Process()
    p.executableURL = URL(fileURLWithPath: "/usr/bin/say")
    p.arguments = [text]
    try? p.run()
}

// ---- TCP client (BSD sockets), with retry while the hub binds ----
func sendAll(_ s: Int32, _ data: [UInt8]) -> Bool {
    var off = 0
    while off < data.count {
        let n = data.withUnsafeBytes { p -> Int in
            send(s, p.baseAddress! + off, data.count - off, 0) }
        if n <= 0 { return false }
        off += n
    }
    return true
}
func sendLog(_ s: Int32, _ id: UInt32, _ text: String) {
    let body = Array(text.utf8)
    var m: [UInt8] = [0x4D, 0x56, 0x4C, 0x47]
    func le32(_ v: UInt32) {
        m.append(UInt8(v & 255)); m.append(UInt8((v >> 8) & 255))
        m.append(UInt8((v >> 16) & 255)); m.append(UInt8((v >> 24) & 255)) }
    le32(id); le32(UInt32(body.count)); m.append(contentsOf: body)
    _ = sendAll(s, m)
}
func readAcks(_ g: Grabber, _ s: Int32) {
    var buf = [UInt8](repeating: 0, count: 4096)
    let n = buf.withUnsafeMutableBytes { recv(s, $0.baseAddress, 4096, MSG_DONTWAIT) }
    if n < 12 { return }
    var i = 0
    while i + 12 <= n {
        if buf[i]==0x4D && buf[i+1]==0x56 && buf[i+2]==0x41 && buf[i+3]==0x4B {
            func u32(_ o: Int) -> UInt32 {
                UInt32(buf[i+o]) | (UInt32(buf[i+o+1])<<8)
                | (UInt32(buf[i+o+2])<<16) | (UInt32(buf[i+o+3])<<24) }
            g.hubRx = u32(4); g.hubDecoded = u32(8)
            g.lastAck = ProcessInfo.processInfo.systemUptime
            i += 12
        } else { i += 1 }
    }
}
func connectLocal(_ port: UInt16) -> Int32 {
    var addr = sockaddr_in()
    addr.sin_family = sa_family_t(AF_INET)
    addr.sin_port = port.bigEndian
    inet_pton(AF_INET, "127.0.0.1", &addr.sin_addr)
    for _ in 0..<50 {                        // ~5 s of retries
        let s = socket(AF_INET, SOCK_STREAM, 0)
        let ok = withUnsafePointer(to: &addr) {
            $0.withMemoryRebound(to: sockaddr.self, capacity: 1) {
                connect(s, $0, socklen_t(MemoryLayout<sockaddr_in>.size)) == 0 } }
        if ok {
            var one: Int32 = 1
            setsockopt(s, IPPROTO_TCP, TCP_NODELAY, &one,
                       socklen_t(MemoryLayout<Int32>.size))
            return s
        }
        close(s); usleep(100_000)
    }
    return -1
}
let sock = connectLocal(port)
guard sock >= 0 else {
    FileHandle.standardError.write(Data("camera 1 could not reach the local hub\n".utf8))
    speak("camera 1 could not reach the local hub")
    hub.interrupt(); exit(1)
}

// ---- camera 1 capture + preview (same as stream_cam) ----
final class Grabber: NSObject, AVCaptureVideoDataOutputSampleBufferDelegate {
    var seq: UInt64 = 0
    var lastSent = 0.0
    var dims = "?"
    var hubRx: UInt32 = 0
    var hubDecoded: UInt32 = 0
    var lastAck = 0.0
    let minGap: Double, sock: Int32
    init(minGap: Double, sock: Int32) { self.minGap = minGap; self.sock = sock }
    func captureOutput(_ o: AVCaptureOutput, didOutput sb: CMSampleBuffer,
                       from c: AVCaptureConnection) {
        let now = ProcessInfo.processInfo.systemUptime
        if now - lastSent < minGap { return }
        lastSent = now
        guard let pb = CMSampleBufferGetImageBuffer(sb) else { return }
        CVPixelBufferLockBaseAddress(pb, .readOnly)
        defer { CVPixelBufferUnlockBaseAddress(pb, .readOnly) }
        let w = CVPixelBufferGetWidth(pb), h = CVPixelBufferGetHeight(pb)
        let stride = CVPixelBufferGetBytesPerRow(pb)
        guard let base = CVPixelBufferGetBaseAddress(pb) else { return }
        var gray = [UInt8](repeating: 0, count: w * h)
        let src = base.assumingMemoryBound(to: UInt8.self)
        for y in 0..<h {
            let row = src + y * stride
            for x in 0..<w {
                let b = Int(row[4*x]), g = Int(row[4*x+1]), r = Int(row[4*x+2])
                gray[y*w+x] = UInt8((77*r + 150*g + 29*b) >> 8)
            }
        }
        var msg = [UInt8](); msg.reserveCapacity(32 + w*h)
        msg.append(contentsOf: [0x4D,0x56,0x46,0x52])
        func le32(_ v: UInt32) {
            msg.append(UInt8(v & 255)); msg.append(UInt8((v>>8) & 255))
            msg.append(UInt8((v>>16) & 255)); msg.append(UInt8((v>>24) & 255)) }
        le32(camid); le32(UInt32(w)); le32(UInt32(h))
        var s = seq; for _ in 0..<8 { msg.append(UInt8(s & 255)); s >>= 8 }
        var t = now.bitPattern; for _ in 0..<8 { msg.append(UInt8(t & 255)); t >>= 8 }
        msg.append(contentsOf: gray)
        if !sendAll(sock, msg) {
            speak("camera 1 lost the local hub")
            exit(1)
        }
        seq += 1; dims = "\(w)x\(h)"
        if seq == 1 { sendLog(sock, camid, "camera 1 (hub-local) connected, \(w)x\(h)") }
        if seq % 50 == 0 { sendLog(sock, camid, "sent \(seq) frames") }
    }
}
let session = AVCaptureSession()
session.sessionPreset = .hd1280x720
guard let dev = AVCaptureDevice.default(for: .video),
      let input = try? AVCaptureDeviceInput(device: dev) else {
    FileHandle.standardError.write(Data("no camera (check permission)\n".utf8))
    speak("no camera on the hub laptop, check camera permission")
    hub.interrupt(); exit(1)
}
session.addInput(input)
let out = AVCaptureVideoDataOutput()
out.videoSettings =
    [kCVPixelBufferPixelFormatTypeKey as String: kCVPixelFormatType_32BGRA]
let grabber = Grabber(minGap: 1.0/fps, sock: sock)
out.setSampleBufferDelegate(grabber, queue: DispatchQueue(label: "grab"))
session.addOutput(out)
session.startRunning()

// on quit, stop the hub child cleanly (so it writes summary.txt)
final class AppDelegate: NSObject, NSApplicationDelegate {
    func applicationShouldTerminateAfterLastWindowClosed(_ s: NSApplication) -> Bool { true }
    func applicationWillTerminate(_ n: Notification) {
        hub.interrupt(); hub.waitUntilExit()
    }
}
let app = NSApplication.shared
let delegate = AppDelegate()
app.delegate = delegate
app.setActivationPolicy(.regular)
let win = NSWindow(contentRect: NSRect(x: 80, y: 80, width: 720, height: 460),
                   styleMask: [.titled, .closable, .resizable, .miniaturizable],
                   backing: .buffered, defer: false)
win.title = "multiview HUB + camera 1  (hub log in the Terminal)"
let preview = AVCaptureVideoPreviewLayer(session: session)
preview.videoGravity = .resizeAspect
let root = NSView(frame: win.contentRect(forFrameRect: win.frame))
root.wantsLayer = true
root.layer = preview
root.autoresizingMask = [.width, .height]
let banner = NSTextField(labelWithString: "starting hub…")
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
var saidPermission = false
let tStart = ProcessInfo.processInfo.systemUptime
Timer.scheduledTimer(withTimeInterval: 0.5, repeats: true) { _ in
    readAcks(grabber, sock)
    if grabber.seq == 0 {
        banner.textColor = .systemYellow
        banner.stringValue = "waiting for camera… (grant camera permission if prompted)"
        if !saidPermission
           && ProcessInfo.processInfo.systemUptime - tStart > 5 {
            saidPermission = true
            speak("hub laptop is waiting for camera permission")
        }
    } else {
        banner.textColor = .systemGreen
        banner.stringValue = "● LIVE  camera 1  \(grabber.dims)  ·  hub got "
            + "\(grabber.hubRx), decoded \(grabber.hubDecoded)  ·  hub IP + "
            + "[cal]/[ext] are in the Terminal"
    }
}
app.run()
