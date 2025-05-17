# q-rs

 * Quantum Computation Simulator for Rust

## Example

### Shor's factoring algorithm

```rust
let n: u32 = 15;
let a = number::coprime(n);
let t: u32 = 3; // precision bits

loop {
    let mut qsim = quantum::Q::new();
    let r0 = qsim.zeros(t);
    let r1 = qsim.zero_log2(n);

    qsim.x(&[r1[r1.len() - 1]]);
    qsim.h(&r0);
    qsim.cmodexp2(a, n, &r0, &r1);
    qsim.iqft(&r0);

    let mut found = false;
    for state in qsim.state().iter() {
        let m0 = state.to_binary_chars(&r0);

        let (s, r, ok) = number::find_order(a, n, &m0);
        if !ok || number::is_odd(r) {
             continue;
        }

        let p0 = number::gcd(a.pow(r / 2) - 1, n);
        let p1 = number::gcd(a.pow(r / 2) + 1, n);
        if number::is_trivial(n, &[p0, p1]) {
            continue;
        }

        println!("{}; s/r={:>2}/{:>2}; p={}, q={}", state, s, r, p0, p1);
        return;
    }
}
```