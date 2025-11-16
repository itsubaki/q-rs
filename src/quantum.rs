use num::Complex;
use num::Zero;

pub type Complex64 = Complex<f64>;

pub type Qubit = Vec<Complex64>;

pub type Gate = Vec<Vec<Complex64>>;

pub type BinaryChars = Vec<char>;

pub struct State {
    number_of_qubits: u32,
    pub index: usize,
    pub amp: Complex64,
    pub prob: f64,
}

impl State {
    pub fn to_binary_chars(&self, qb: &[u32]) -> BinaryChars {
        let v = to_binary_chars(self.index, self.number_of_qubits as usize);

        let mut bin = vec![];
        for i in qb {
            bin.push(v[*i as usize]);
        }

        bin
    }
}

impl std::fmt::Display for State {
    fn fmt(&self, dest: &mut std::fmt::Formatter) -> std::fmt::Result {
        let bits: String = format!("{:>0n$b}", self.index, n = self.number_of_qubits as usize);
        write!(
            dest,
            "[{}]({:>+.4} {:>+.4}): {:>.4}",
            bits, self.amp.re, self.amp.im, self.prob,
        )
    }
}

pub struct Q {
    qb: Qubit,
}

impl Q {
    pub fn new() -> Q {
        Q { qb: vec![] }
    }

    pub fn add(&mut self, qb: Qubit) -> u32 {
        if self.qb.is_empty() {
            self.qb = qb;
            return 0;
        }

        self.tensor(qb);
        self.number_of_qubits() - 1
    }

    pub fn zero(&mut self) -> u32 {
        self.add(vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)])
    }

    pub fn zeros(&mut self, n: u32) -> Vec<u32> {
        let mut list = vec![];

        for _ in 0..n {
            list.push(self.zero());
        }

        list
    }

    pub fn zero_log2(&mut self, n: u32) -> Vec<u32> {
        let log2n = ((n as f64).log2() as u32) + 1;
        self.zeros(log2n)
    }

    fn tensor(&mut self, qb: Qubit) {
        let mut v: Qubit = vec![];

        for x in &self.qb {
            for y in &qb {
                v.push(x * y);
            }
        }

        self.qb = v
    }

    pub fn number_of_qubits(&self) -> u32 {
        (self.qb.len() as f64).log2() as u32
    }

    pub fn x(&mut self, qb: &[u32]) {
        self.apply_with(x(), qb)
    }

    pub fn h(&mut self, qb: &[u32]) {
        self.apply_with(h(), qb)
    }

    fn apply_with(&mut self, g: Gate, qb: &[u32]) {
        let list: Vec<Gate> = gate_list(self.number_of_qubits(), g, qb);
        let g: Gate = tensor_with(&list);

        self.apply(g)
    }

    pub fn apply(&mut self, g: Gate) {
        let mut v: Qubit = vec![];

        for item in &g {
            let mut e = Complex::new(0.0, 0.0);

            for (j, _) in item.iter().enumerate() {
                e += item[j] * self.qb[j];
            }

            v.push(e);
        }

        self.qb = v
    }

    pub fn iqft(&mut self, qb: &[u32]) {
        let len = qb.len();

        for i in (0..len).rev() {
            let mut k = (len - i) as i32;

            for j in ((i + 1)..len).rev() {
                let theta = -2.0 * std::f64::consts::PI / (2.0_f64.powf(k as f64));
                self.cr(theta, qb[j], qb[i]);
                k -= 1;
            }

            self.h(&[qb[i]]);
        }
    }

    pub fn cr(&mut self, theta: f64, control: u32, target: u32) {
        let n = self.number_of_qubits();
        let g: Gate = cr(theta, n, control, target);
        self.apply(g)
    }

    pub fn state(&self) -> Vec<State> {
        let mut list = vec![];
        let nob = self.number_of_qubits();

        for i in 0..self.qb.len() {
            let r = round(self.qb[i]);
            if r.is_zero() {
                continue;
            }

            list.push(State {
                number_of_qubits: nob,
                index: i,
                amp: r,
                prob: r.norm().powf(2.0),
            });
        }

        list
    }
}

impl Default for Q {
    fn default() -> Self {
        Self::new()
    }
}

fn gate_list(nob: u32, g: Gate, qb: &[u32]) -> Vec<Gate> {
    let mut list: Vec<Gate> = vec![];

    for i in 0..nob {
        let mut found = false;

        for j in qb {
            if i == *j {
                found = true;
                break;
            }
        }

        if found {
            list.push(g.to_vec());
            continue;
        }

        list.push(id(1));
    }

    list
}

fn tensor_with(list: &[Gate]) -> Gate {
    let mut g: Gate = list[0].to_vec();

    for i in list.iter().skip(1) {
        g = tensor(g, i.to_vec());
    }

    g
}

fn tensor(m: Gate, n: Gate) -> Gate {
    let mut g: Gate = vec![];

    for (i, _) in m.iter().enumerate() {
        for (k, _) in n.iter().enumerate() {
            let mut v = vec![];

            for (j, _) in m[i].iter().enumerate() {
                for (l, _) in n[k].iter().enumerate() {
                    v.push(m[i][j] * n[k][l]);
                }
            }

            g.push(v);
        }
    }

    g
}

pub fn x() -> Gate {
    vec![
        vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
        vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
    ]
}

pub fn h() -> Gate {
    let e = Complex::new(1.0 / std::f64::consts::SQRT_2, 0.0);
    vec![vec![e, e], vec![e, -1.0 * e]]
}

pub fn id(nob: u32) -> Gate {
    let s = 1 << nob;
    let mut g = vec![vec![Complex::new(0.0, 0.0); s]; s];

    for (i, row) in g.iter_mut().enumerate() {
        row[i] = Complex::new(1.0, 0.0);
    }

    g
}

pub fn cr(theta: f64, nob: u32, control: u32, target: u32) -> Gate {
    let mut g: Gate = id(nob);
    let e = Complex::new(0.0, theta).exp();
    let mask = (1 << (nob - 1 - control)) as usize;

    for (i, v) in g.iter_mut().enumerate() {
        if (i & mask) == mask && (i & (1 << (nob - 1 - target))) != 0 {
            v[i] = e * v[i];
        }
    }

    g
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

fn to_binary_chars(i: usize, nob: usize) -> BinaryChars {
    format!("{:>0n$b}", i, n = nob).chars().collect()
}

#[test]
fn test_to_binary_chars() {
    assert_eq!(to_binary_chars(3, 5), vec!['0', '0', '0', '1', '1']);
    assert_eq!(to_binary_chars(7, 5), vec!['0', '0', '1', '1', '1']);
    assert_eq!(to_binary_chars(15, 5), vec!['0', '1', '1', '1', '1']);
    assert_eq!(to_binary_chars(31, 5), vec!['1', '1', '1', '1', '1']);
}
