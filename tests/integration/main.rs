//! Every integration test, in one test binary.
//!
//! Cargo runs separate test binaries one after another, and each runs its own
//! tests in parallel — so with six binaries the suite's wall time was the SUM of
//! six slowest-test tails, and the machine sat at 3.6 of 10 cores. One binary
//! puts every test in a single pool, so wall time is bounded by the slowest
//! single test or by total work over the core count, whichever is larger. It
//! also links once instead of six times.
//!
//! Run one area with a module filter, e.g.
//! `cargo test --release --test integration end_to_end::`.
//!
//! **Module names set the schedule.** libtest starts tests in name order and has
//! no stable way to change that, so a slow test that sorts last starts last and
//! the run ends on one or two busy cores. The whole-pipeline compare tests
//! (formerly `stage4_integration`) hold most of the slowest single tests, so
//! they live in `end_to_end`, which sorts before `stage*`: simulated from
//! measured durations, 258 s against 285 s for the old order and 256 s for an
//! ideal longest-first schedule. Keep slow end-to-end tests in that module.

mod end_to_end;
mod stage1_integration;
mod stage2_competition;
mod stage2_integration;
mod stage2_segmentation;
mod stage3_integration;
