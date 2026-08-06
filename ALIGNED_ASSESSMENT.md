# Assessment: the `aligned` project, and what it implies for rig solving

Date of assessment: **2026-08-06**. This note records a design discussion so it
does not have to be repeated. It is a decision record, not a work log; see
`CONVERSATION_LOG.md` for session history.

The short version: `aligned` is **not** a drop-in for anything here, but the
*architecture* it embodies is the right answer to a real defect in our
correspondence logic, and we intend to refactor along those lines.

---

## 1. What `aligned` is

`aligned` is a sibling project at `../aligned`, independent of this repo. Its
`NOTES.md:17-33` states its purpose: a portable C99 library for **temporally
aligning camera frames captured on two different Macs**. Two laptops joined by
USB-C, `BA` (conductor, machine MA) and `BB` (follower, MB); MB streams frames
to MA; MA pairs them by an estimate of **when photons hit each sensor**, not
when the frame was delivered.

Its pieces: an NTP-style clock model, an optimal non-crossing dynamic-programming
(DP) matcher, a preallocated slab frame store, a transport vtable, and a reader
API for downstream apps.

**State as of 2026-08-06:** step 1 of a 9-step build plan. Zero `.c` files. The
only artifact is `include/aligned/aligned.h` (123 lines, covering status, time,
and clock only). *This section will go stale — re-check before relying on it.*

`aligned`'s own `NOTES.md:32-33` declares it independent of every sibling repo,
in both directions. Any use of it here would need that constraint revisited with
its owner.

## 2. Why it is not a drop-in

Recorded so this is not re-argued:

1. **It contains no geometry.** "Rig solve" here means recovering
   `x_a = R x_b + t` (`tools/hub_solve.c`). `aligned` has no rotation, pose, or
   SE(3) concept anywhere. It solves *temporal* alignment. The overlap is with
   `tools/hub_clock.c` and the pairing loop, never with the solver.
2. **Where it overlaps, ours is ahead.** See §5.
3. **Shape mismatch.** `aligned` is two hosts with MA as sole conductor driving
   capture trials; `livehub.c` already owns the socket, the session state
   machine, and the frozen MVPB/MVTS wire format across up to 8 cameras.

## 3. The real finding: our correspondence gate is not honest

This is the substantive outcome of the discussion, and it is a defect in *this*
repo, independent of whether `aligned` is ever used.

The pairing loop at `tools/livehub.c:495-533` decides which anchor observations
from two cameras describe the same physical instant. In its relaxed path it
accepts:

- **up to 15 seconds** of clock skew between the two observations
  (`livehub.c:511`), and
- **up to 120 counter ticks** of Gray-code difference (`livehub.c:521`),

both unlocked by the `dwell` flag, which means only "the display looked static
compared to the previous anchor." Additionally, the coarse tier compares
counters mod 256 (`livehub.c:516-519`), so two genuinely different display
states 256 apart are indistinguishable and are treated as a match.

So in relaxed mode we are **not** pairing simultaneous observations. We pair
things up to 15 s apart and assert that nothing moved in between.

**Why the relaxation exists** — and this matters, because it was not
carelessness. Anchors arrive every 7–15 s per camera (`livehub.c:509-510`). Two
cameras with anchor spacing that long and no phase relationship will almost
never produce a strictly simultaneous pair. Given the architecture, relaxing was
the only way to get observations at all.

**Why it is harmful.** `hub_solve_pair()` treats each `hub_obs` as an
independent measurement of one fixed rigid transform. If the rig moved during
those 15 seconds, the observation is not noisy — it is measuring a *different
quantity*. That is exactly the failure mode `tools/hub_solve.c:5-11` describes
and defends against: coherent outliers that all pull the same way. The seeded
Huber-IRLS with two-pass hard trim is sophisticated machinery built to survive a
problem manufactured upstream by the gate. That coupling — a permissive gate
paid for by a robust estimator — is the kind of long-range interaction that is
hard to debug, and it is the reason for the refactor below.

**The Gray-code counter itself is not the problem.** It is a uniquely-decodable
timestamp physically embedded in the observed scene and seen by both cameras
through the same channel — a clapperboard, and a better timebase than two
crystal oscillators talking over TCP. Keep it. The defect is that `dwell` was
allowed to override it whenever the honest answer was inconvenient.

## 4. Decision: hard separation between pair collector and pair consumer

Adopt `aligned`'s **contract and discipline**, not its code.

| | |
|---|---|
| **Below the line** (collector) | Measured time only. No content reasoning, no dwell exceptions, no widening the gate because the count came back low. Emits pairs + Δt + an uncertainty. |
| **Above the line** (consumer) | Gray-code decode, PnP, `hub_obs`, `hub_solve`. |

Rationale: compartmentalize the complexity budget. The rig solver should consume
matched pairs without the counter/tier/dwell details polluting its concept
space, and the collector should be forced to report sensed reality as measured
rather than as massaged.

This also resolves the counter's double duty. Today it serves as both the
correspondence signal *and* the thing being observed. Under the split, the first
role disappears entirely (measured time does that job) and the second stays above
the line, where it belongs.

Two consequences worth having:

- `hub_solve`'s robustness machinery becomes **insurance rather than
  load-bearing**. Do not delete it; stop depending on it to paper over bad
  correspondence.
- **Replay stops producing misleading numbers.** `replaycam` time-scrambles
  counter matches, which is why `CONVERSATION_LOG.md` must caveat that replay
  rig estimates are not accuracy benchmarks. Under a temporal contract, replay
  either supplies honest capture times or the collector reports zero pairs.

## 5. Decision: keep our clock, make the clock method pluggable

"Our clock" is `tools/hub_clock.c` — strictly not a clock but a per-camera model
`hub_time = phi + rate * t_cam` (`hub_clock.c:10`), fit from MVPB/MVTS probes at
~1 Hz, mapping each camera's capture timestamps into the hub's `CLOCK_MONOTONIC`
domain. Fit is Theil-Sen slope + median intercept, then a MAD-gated least-squares
refinement on inliers (`hub_clock.c:23`).

Prefer it over `aligned`'s planned min-RTT + weighted least squares because its
gates (`hub_clock.c:16-22`) each exist because live Wi-Fi produced the specific
misbehavior guarded against: RTT over 25 ms or 3× the ring median excluded;
fewer than 5 accepted samples means no fit; fitted `|rate-1|` over 1% rejected as
garbage; `t_cam` moving backwards resets that camera's history. Theil-Sen also
has a far better breakdown point than WLS. Under a hard temporal boundary the
clock becomes load-bearing, so this matters more, not less.

**Caveats on trusting it.** There is an identity fallback: a camera that never
answers MVPB, or whose fit fails a gate, maps `t` to itself, i.e. reverts to hub
receive time. `livehub.c:493-494` checks `hub_clock_err() >= 0` for both cameras
first, so the pairing loop knows the difference — but "synced" is a per-camera
runtime property, not a guarantee. Also, the 0.16–0.32 ms figure in
`CONVERSATION_LOG.md` is from synthetic tests in `test_clock_sync`; the live
claim in the source header is the softer ~1–2 ms.

`aligned`'s clock API is already nearly the right shape — `..._model_fit()` plus
the two domain conversions (`aligned.h:66-117`). Making it a vtable is consistent
with the transport and frame-source vtables it already has.

**One naming constraint** (cheap now, painful later): the conversion should be
`to_reference(host, t)`, not `to_ma`/`to_mb`. Every backend implements against
that signature, so it is the one piece that is expensive to change later.

## 6. The key insight: pair-then-decode, not decode-then-pair

This defeats the obvious objection that strict gating will collapse the pair
count.

- **Today (decode-then-pair):** decode opportunistically, then hope two decoded
  anchors happen to be simultaneous. With 7–15 s anchor spacing and no phase
  relationship, a strict window almost never fires. Hence §3's relaxation.
- **Proposed (pair-then-decode):** capture everything at video rate, store it,
  match temporally on timestamps alone (cheap), and spend the expensive decode
  only on frames that already have a simultaneous partner.

Two free-running 30 fps sensors with arbitrary relative phase have a
nearest-partner |Δt| averaging about a quarter of the frame interval — roughly
8 ms, worst case ~17 ms. Essentially *every* frame has a genuinely simultaneous
partner. Pair count stops being a matter of luck and becomes a matter of capture
rate and decode budget, with the entire decode budget spent on frames that can
actually yield a rig observation. The reader is now ~17× faster and the pipeline
is already feed-limited rather than decode-limited, so this is the right way
round.

`aligned`'s phase nudging (`NOTES.md:48`, restarting MB's session by the measured
residual phase) attacks the residual directly, improving pair quality by
experiment across trials — the honest version of what dwell relaxation faked.

## 7. Decision: get N=2 right before generalizing

Build and validate the two-camera case against real hardware first; pay later to
raise N. Do not generalize to N collectors now.

Rationale: the rig consumer is *already* pairwise — `hub_solve_pair()` solves one
`(a,b)`, and `livehub.c` carries `rig_a`/`rig_b` and `[MAXCAMS][MAXCAMS]` pair
state. The live 0.319 m baseline is a two-camera number. An N=2 collector matches
the unit the consumer actually works in, and a general design never validated
against reality is a worse foundation than a specific one that has been. The
`to_reference` naming in §5 is the only concession made to future N, and it costs
nothing today.

## 8. Open question / next experiment

**Measure the frame-arrival Δt distribution between two cameras at full capture
rate, before any decode.**

This is the number the whole design rests on, it requires no refactor to obtain,
and it decides whether §6's phase argument holds in practice. If it comes back
where the arithmetic says it should, the pair-count concern is dead and the
refactor is justified on measurement rather than argument. If it does not, the
real problem is capture-side and we want to know that first.

Note this supersedes an earlier, worse experiment (tighten the existing gate and
observe): that measures only the sparse decode-then-pair regime we intend to
leave.

## 9. Considered and rejected

- **Vendoring or depending on `aligned` as a library.** Rejected: no
  implementation exists yet, shape mismatch (§2), and it would add a second
  clock model, transport, and frame store — spending complexity budget rather
  than saving it.
- **Layering over a *frame-pair* boundary while keeping the counter as the
  correspondence signal.** Rejected: keeps the massage in place behind a new
  abstraction, giving both the abstraction and the thing it was meant to hide.
- **Adopting `aligned`'s clock model.** Rejected in favour of `hub_clock`
  (§5), but its *pluggability* is adopted.
