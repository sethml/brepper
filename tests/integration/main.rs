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
//! `cargo test --release --test integration stage4_integration::`.

mod stage1_integration;
mod stage2_competition;
mod stage2_integration;
mod stage2_segmentation;
mod stage3_integration;
mod stage4_integration;
