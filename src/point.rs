// use core::ops::{Add, AddAssign, Mul, Neg, Sub};
// use num_traits::{One, Zero};

use crate::field::FieldElement;

// 定义 Point 结构体
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: FieldElement,
    pub y: FieldElement,
}

impl Point {
    // 创建新的 Point
    pub fn new(x: FieldElement, y: FieldElement) -> Self {
        debug_assert_eq!(x.square() + y.square(), FieldElement::one());
        Point { x, y }
    }
    pub fn zero() -> Self {
        //生成零点（1,0)
        Self::new(FieldElement::one(), FieldElement::zero())
    }
    pub fn double(self) -> Self {
        Self::new(
            self.x.square().double() - FieldElement::one(),
            self.x.double() * self.y,
        )
    }
    pub fn add(self, other: Self) -> Self {
        Self::new(
            self.x * other.x - self.y * other.y,
            self.x * other.y + other.x * self.y,
        )
    }

    pub fn npower(&self, mut n: i64) -> Self {
        let mut result = Self::zero();
        let mut a = self.clone();
        while n > 0 {
            if n % 2 == 1 {
                result = result.add(a)
            }
            a = a.double();
            n = n / 2;
        }
        result
    }
    pub fn get_modulus() -> i64 {
        let zero = Point::zero();
        let modulus = zero.x.get_modulus().clone();
        modulus
    }
    //
    pub fn if_primitive_generator(&self) -> bool {
        let p = Point::get_modulus();
        let a = self.npower(p + 1) == Self::zero() && self.npower((p + 1) / 2) != Self::zero();
        a
    }

    //get the primitive_generator of C(Fp):find point and check if it is primitive_generator.
    pub fn get_primitive_generator() -> Self {
        let mut g = Self::new(FieldElement::one(), FieldElement::zero());
        let p = g.x.get_modulus().clone();
        'outer: for x in 1..p - 1 {
            'inner: for y in 1..p - 1 {
                if (x * x + y * y) % p == 1 {
                    g = Self::new(FieldElement::new(x), FieldElement::new(y));
                    if g.if_primitive_generator() {
                        break 'outer;
                    }
                }
            }
        }
        g
    }
    //get the generator of Gn(n actually is logn)
    pub fn get_nth_generator(logn: u64) -> Self {
        let g = Point::get_primitive_generator();
        let p = g.x.get_modulus().clone();
        debug_assert!(
            1 << logn <= p + 1,
            "Debug assertion failed: n should be less than p"
        );
        let n = (p + 1) / (1 << logn);
        g.npower(n)
    }
    pub fn generate_Gn_byG(G: Point, logn: u64) -> Vec<Point> {
        let mut result: Vec<Point> = Vec::new();
        for i in 0..(1 << logn) {
            let point = G.npower(i);
            result.push(point);
        }
        result
    }

    pub fn spcoset(logn: u64) -> Vec<Point> {
        let Gn = Point::get_nth_generator(logn); //the generator of Gn
        let Q = Point::get_nth_generator(logn + 1); //the shift，and the order should be 2^(logn+1)
        let mut result: Vec<Point> = Vec::new();
        for i in 0..(1 << logn) {
            let Gnpoint = Q.add(Gn.npower(i));
            result.push(Gnpoint);
        }
        result
    }
    pub fn spcoset_by_Q(Q: Point, logn: u64) -> Vec<Point> {
        let Gn = Point::get_nth_generator(logn);
        let mut result: Vec<Point> = Vec::new();
        for i in 0..(1 << logn) {
            let Gnpoint = Q.add(Gn.npower(i));
            result.push(Gnpoint);
        }
        result
    }

    pub fn v_n(mut self, log_n: usize) -> FieldElement {
        for _ in 0..(log_n - 1) {
            self.x = self.x.square().double() - FieldElement::one(); //即：计算vn(x)=2*vn-1(x)^2-1。这里迭代log_n次，推测应是N=2^(log_n)
        }
        self.x
    }
}

// 示例用法
#[test]
fn main() {
    // Generate Gn;
    let g_5 = Point::new(FieldElement::new(10), FieldElement::new(5));
    let g_4 = Point::new(FieldElement::new(13), FieldElement::new(7));
    let g_3 = Point::new(FieldElement::new(27), FieldElement::new(27));
    let cfp = Point::generate_Gn_byG(g_5, 5);
    let G4 = Point::generate_Gn_byG(g_4, 4);
    let G3 = Point::generate_Gn_byG(g_3, 3);
    let G3_prime = Point::generate_Gn_byG(g_4, 3);
    for (index, value) in G3.iter().enumerate() {
        println!("Index:{:?}, Value: {:?}", index, value);
    }
    println!("get subgroup Gn!");
    println!("------------------------------------");

    // Generate standard position twin-cosets;
    let Q = Point::new(FieldElement::new(24), FieldElement::new(13));
    let spcoset = Point::spcoset_by_Q(Q, 3);
    for (index, value) in spcoset.iter().enumerate() {
        println!("Index:{:?}, Value: {:?}", index, value);
    }
    println!("get standard position cosets D!");
    println!("------------------------------------");
}
