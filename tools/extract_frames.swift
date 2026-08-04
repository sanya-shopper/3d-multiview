// multiview frame extractor: HEVC/H.264 video -> 8-bit gray PGM frames.
// Uses only macOS system frameworks (AVFoundation/CoreGraphics), keeping
// the no-external-dependency rule. Usage:
//   swiftc -O tools/extract_frames.swift -o /tmp/extract_frames
//   /tmp/extract_frames <in.mov> <outdir> <interval-seconds> [start-s end-s]
// The optional window limits extraction to [start, end] seconds -- used
// to mine dense frames around instants chosen via the display counter.
import AVFoundation
import CoreGraphics
import Foundation

let a = CommandLine.arguments
guard a.count == 4 || a.count == 6, let interval = Double(a[3]) else {
    print("usage: extract_frames <in.mov> <outdir> <interval-s> [start-s end-s]")
    exit(1)
}
let asset = AVURLAsset(url: URL(fileURLWithPath: a[1]))
let gen = AVAssetImageGenerator(asset: asset)
gen.requestedTimeToleranceBefore = .zero
gen.requestedTimeToleranceAfter = .zero
gen.appliesPreferredTrackTransform = true
var dur = CMTimeGetSeconds(asset.duration)
var t = 0.25, idx = 0
if a.count == 6, let ws = Double(a[4]), let we = Double(a[5]) {
    t = ws
    if we < dur { dur = we }
}
while t < dur {
    let time = CMTime(seconds: t, preferredTimescale: 600)
    if let cg = try? gen.copyCGImage(at: time, actualTime: nil) {
        let w = cg.width, h = cg.height
        var buf = [UInt8](repeating: 0, count: w * h)
        let cs = CGColorSpaceCreateDeviceGray()
        buf.withUnsafeMutableBytes { p in
            if let ctx = CGContext(data: p.baseAddress, width: w,
                                   height: h, bitsPerComponent: 8,
                                   bytesPerRow: w, space: cs,
                bitmapInfo: CGImageAlphaInfo.none.rawValue) {
                ctx.draw(cg, in: CGRect(x: 0, y: 0, width: w, height: h))
            }
        }
        var data = "P5\n\(w) \(h)\n255\n".data(using: .ascii)!
        data.append(contentsOf: buf)
        try? data.write(to: URL(fileURLWithPath:
            String(format: "%@/f%04d.pgm", a[2], idx)))
        idx += 1
    }
    t += interval
}
print("extracted \(idx) frames (\(Int(dur))s) from \(a[1])")
