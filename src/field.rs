use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldElement(i64); // 元组结构体

impl FieldElement {
    const P: i64 = 31; // 模数

    pub fn new(value: i64) -> Self {
        assert!(
            Self::P > 0 && Self::P % 4 == 3,
            "p must be a positive integer and satisfy p % 4 == 3"
        );
        FieldElement((value % Self::P + Self::P) % Self::P)
    }
    pub fn get_modulus(&self) -> i64 {
        FieldElement::P
    }

    pub fn inverse(&self) -> Option<Self> {
        if self.0 == 0 {
            return None;
        }
        let (mut a, mut b) = (self.0, Self::P);
        let (mut x, mut y) = (1, 0);
        while b != 0 {
            let q = a / b;
            (a, b) = (b, a - q * b);
            (x, y) = (y, x - q * y);
        }

        if a == 1 {
            Some(FieldElement::new(x))
        } else {
            None
        }
    }

    pub fn one() -> Self {
        FieldElement::new(1)
    }

    pub fn zero() -> Self {
        FieldElement::new(0)
    }

    pub fn square(self) -> Self {
        self.clone() * self.clone()
    }
    pub fn double(&self) -> Self {
        self.clone() + self.clone()
    }
}

// 实现 From trait 以允许直接从 i64 转换
impl From<i64> for FieldElement {
    fn from(value: i64) -> Self {
        FieldElement::new(value)
    }
}

// 直接实现运算符
impl Add for FieldElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        FieldElement::new(self.0 + rhs.0)
    }
}

impl Sub for FieldElement {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        FieldElement::new(self.0 - rhs.0)
    }
}

impl Mul for FieldElement {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        FieldElement::new(self.0 * rhs.0)
    }
}

impl Div for FieldElement {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inverse().unwrap()
    }
}

impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// 通过 Deref Trait 实现直接使用 FieldElement(3)
impl std::ops::Deref for FieldElement {
    type Target = i64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[test]
fn test() {
    let a = FieldElement(9); // 仍然可以使用这种方式
    let b = FieldElement(5);
    let c = FieldElement::one();
    let d = FieldElement::zero();
    // 现在可以直接用 FieldElement(3) 来创建实例
    let e = b.square();

    let g = c.square();

    println!("e = {e}");
    println!("a * c = {}", a * c);
}
