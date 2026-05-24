use num::Complex;
use num::Zero;

pub type Complex64 = Complex<f64>;

pub type Qubit = Vec<Complex64>;

pub type Gate = Vec<Vec<Complex64>>;

pub type BinaryChars = Vec<char>;

pub struct State {
    number_of_qubits: usize,
    pub index: usize,
    pub amp: Complex64,
    pub prob: f64,
}

impl State {
    pub fn to_binary_chars(&self, qb: &[usize]) -> BinaryChars {
        let v = to_binary_chars(self.index, self.number_of_qubits);

        let mut bin = vec![];
        for i in qb {
            bin.push(v[*i]);
        }

        bin
    }
}

impl std::fmt::Display for State {
    fn fmt(&self, dest: &mut std::fmt::Formatter) -> std::fmt::Result {
        let bits: String = format!("{:>0n$b}", self.index, n = self.number_of_qubits);
        write!(
            dest,
            "[{}]({:>+.4} {:>+.4}): {:>.4}",
            bits, self.amp.re, self.amp.im, self.prob,
        )
    }
}

pub struct Q {
    qb: Qubit,
    number_of_qubits: usize,
}

impl Default for Q {
    fn default() -> Self {
        Self {
            qb: Qubit::new(),
            number_of_qubits: 1,
        }
    }
}

impl Q {
    pub fn new() -> Self {
        Self {
            qb: vec![],
            number_of_qubits: 0,
        }
    }

    pub fn number_of_qubits(&self) -> usize {
        self.number_of_qubits
    }

    pub fn zero(&mut self) -> usize {
        if self.qb.is_empty() {
            self.qb = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];

            self.number_of_qubits += 1;
            return 0;
        }

        let now = std::mem::take(&mut self.qb);
        self.qb = Vec::with_capacity(now.len() * 2);
        for v in now {
            self.qb.push(v);
            self.qb.push(Complex64::new(0.0, 0.0));
        }

        self.number_of_qubits += 1;
        self.number_of_qubits - 1
    }

    pub fn zeros(&mut self, n: usize) -> Vec<usize> {
        (0..n).map(|_| self.zero()).collect()
    }

    pub fn apply(&mut self, g: Gate) {
        let mut next = Vec::with_capacity(self.qb.len());
        for row in &g {
            let mut acc = Complex64::new(0.0, 0.0);
            for (&gij, &qj) in row.iter().zip(&self.qb) {
                acc += gij * qj;
            }

            next.push(acc);
        }

        self.qb = next;
    }

    fn g(&mut self, target: usize, m00: Complex64, m01: Complex64, m10: Complex64, m11: Complex64) {
        let bit = 1usize << (self.number_of_qubits - 1 - target);

        for i in 0..self.qb.len() {
            if (i & bit) != 0 {
                continue;
            }

            let j = i | bit;
            let a = self.qb[i];
            let b = self.qb[j];
            self.qb[i] = m00 * a + m01 * b;
            self.qb[j] = m10 * a + m11 * b;
        }
    }

    pub fn x(&mut self, qb: &[usize]) {
        for &q in qb {
            self.g(
                q,
                Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0),
            );
        }
    }

    pub fn h(&mut self, qb: &[usize]) {
        let s = 1.0 / std::f64::consts::SQRT_2;
        for &q in qb {
            self.g(
                q,
                Complex64::new(s, 0.0),
                Complex64::new(s, 0.0),
                Complex64::new(s, 0.0),
                Complex64::new(-s, 0.0),
            );
        }
    }

    pub fn cr(&mut self, theta: f64, control: usize, target: usize) {
        let cbit = 1usize << (self.number_of_qubits - 1 - control);
        let tbit = 1usize << (self.number_of_qubits - 1 - target);
        let phase = Complex64::new(0.0, theta).exp();

        for i in 0..self.qb.len() {
            if (i & cbit) != 0 && (i & tbit) != 0 {
                self.qb[i] *= phase;
            }
        }
    }

    pub fn iqft(&mut self, qb: &[usize]) {
        let len = qb.len();
        for i in (0..len).rev() {
            let mut k = (len - i) as i32;
            for j in ((i + 1)..len).rev() {
                let theta = -2.0 * std::f64::consts::PI / (2.0_f64.powi(k));
                self.cr(theta, qb[j], qb[i]);
                k -= 1;
            }

            self.h(&[qb[i]]);
        }
    }

    pub fn state(&self) -> Vec<State> {
        let mut list = vec![];
        for (i, &amp) in self.qb.iter().enumerate() {
            let amp = round(amp);

            if amp.is_zero() {
                continue;
            }

            list.push(State {
                number_of_qubits: self.number_of_qubits,
                index: i,
                amp,
                prob: amp.norm_sqr(),
            });
        }

        list
    }
}

fn round(c: Complex64) -> Complex64 {
    let mut round = c;
    if c.re.abs() < 1e-13 {
        round.re = 0.0;
    }

    if c.im.abs() < 1e-13 {
        round.im = 0.0;
    }

    round
}

fn to_binary_chars(i: usize, num: usize) -> BinaryChars {
    format!("{:>0n$b}", i, n = num).chars().collect()
}

#[test]
fn test_to_binary_chars() {
    assert_eq!(to_binary_chars(3, 5), vec!['0', '0', '0', '1', '1']);
    assert_eq!(to_binary_chars(7, 5), vec!['0', '0', '1', '1', '1']);
    assert_eq!(to_binary_chars(15, 5), vec!['0', '1', '1', '1', '1']);
    assert_eq!(to_binary_chars(31, 5), vec!['1', '1', '1', '1', '1']);
}
