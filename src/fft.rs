// use core::ops::{Add, AddAssign, Mul, Neg, Sub};
// use num_traits::{One, Zero};

use crate::field::FieldElement;
use crate::point::Point;
use crate::utils;

// It's the first layer of ifft, tranfer evaluations to coffecients。 The different between first layer with other layer is that:（1）element of Domain is Point(2)factors is y
fn interpolate(vals: &Vec<FieldElement>, domain: &Vec<Point>) -> Vec<FieldElement> {
    let half_length = domain.len() / 2;
    //first layer：just pick the first half domain points,and project to x
    let half_domain: Vec<FieldElement> = utils::get_halfdomain_of_point(domain);

    //split the vals to two part
    let left: Vec<FieldElement> = vals.get(..half_length).unwrap_or_default().to_vec();
    let mut right: Vec<FieldElement> = vals.get(half_length..).unwrap_or_default().to_vec();
    right.reverse(); //in circle fft, we need let the right part do reverse.
    let mut f0: Vec<FieldElement> = Vec::with_capacity(half_length); //Next left layer
    let mut f1: Vec<FieldElement> = Vec::with_capacity(half_length); //Next right layer

    //factor of first layer: 1/y
    let factor = utils::get_factor_of_point(domain);
    for i in 0..left.len() {
        let l = (left[i] + right[i]) / FieldElement::new(2);
        f0.push(l);
        let r = (left[i] - right[i]) / (factor[i] * FieldElement::new(2));
        f1.push(r);
    }
    //recrusive
    let mut o = vec![FieldElement::zero(); domain.len()];
    let f0_result = ifft(&f0, &half_domain);
    let f1_result = ifft(&f1, &half_domain);
    for i in 0..half_domain.len() {
        o[2 * i] = f0_result[i];
        o[2 * i + 1] = f1_result[i];
    }
    o
}

// It's the regular layer of ifft, tranfer evaluations to coffecients
fn ifft(vals: &Vec<FieldElement>, domain: &Vec<FieldElement>) -> Vec<FieldElement> {
    let n = vals.len();
    if n == 1 {
        return vec![vals[0]];
    }
    let half_length = domain.len() / 2;
    //half domain is the first half part
    let half_domain: Vec<FieldElement> = utils::get_halfdomain_of_vec(domain);

    //split the vals to two part
    let left: Vec<FieldElement> = vals.get(..half_length).unwrap_or_default().to_vec();
    let mut right: Vec<FieldElement> = vals.get(half_length..).unwrap_or_default().to_vec();
    right.reverse(); //in circle fft, we need let the right part do reverse.

    let mut f0: Vec<FieldElement> = Vec::with_capacity(half_length);
    let mut f1: Vec<FieldElement> = Vec::with_capacity(half_length);
    //factor of first layer: 1/x
    let factor = utils::get_factor_of_vec(domain);
    for i in 0..left.len() {
        let l = (left[i] + right[i]) / FieldElement::new(2);
        f0.push(l);
        let r = (left[i] - right[i]) / (factor[i] * FieldElement::new(2));
        f1.push(r);
    }
    //recrusive
    let mut o = vec![FieldElement::zero(); domain.len()];
    let f0_result = ifft(&f0, &half_domain);
    let f1_result = ifft(&f1, &half_domain);
    for i in 0..half_domain.len() {
        o[2 * i] = f0_result[i];
        o[2 * i + 1] = f1_result[i];
    }
    o
}

// It's the first layer of fft, tranfer coffecients to evaluations. The different between first layer with other layer is that:（1）element of Domain is Point(2)factors is y
fn evaluate(coffecients: &Vec<FieldElement>, domain: &Vec<Point>) -> Vec<FieldElement> {
    let half_length = &domain.len() / 2;
    //first layer：just pick the first half domain points,and project to x
    let half_domain: Vec<FieldElement> = utils::get_halfdomain_of_point(domain);

    let evens: Vec<_> = coffecients.iter().step_by(2).cloned().collect();
    let odds: Vec<_> = coffecients.iter().skip(1).step_by(2).cloned().collect();
    let f0 = fft(&evens, &half_domain); //extract the even part,and do fft recursivly
    let f1 = fft(&odds, &half_domain); //extract the odd part,and do fft recursivly

    let mut left: Vec<FieldElement> = Vec::with_capacity(half_length);
    let mut right: Vec<FieldElement> = Vec::with_capacity(half_length);
    let factor = utils::get_factor_of_point(domain);
    for i in 0..factor.len() {
        let l = f0[i] + factor[i] * f1[i];
        left.push(l);
        let r = f0[i] - factor[i] * f1[i];
        right.push(r);
    }
    right.reverse(); //in circle fft, we need let the right part do reverse.

    let mut o = vec![FieldElement::zero(); domain.len()];
    for i in 0..half_domain.len() {
        o[i] = left[i];
        o[i + half_domain.len()] = right[i]; //combine,then it is the evaluations
    }
    return o;
}

// It's the regular layer of fft, tranfer coffecients to evaluations
fn fft(coffecients: &Vec<FieldElement>, domain: &Vec<FieldElement>) -> Vec<FieldElement> {
    if coffecients.len() == 1 {
        return vec![coffecients[0]];
    }
    let half_length = domain.len() / 2;
    let half_domain: Vec<FieldElement> = utils::get_halfdomain_of_vec(domain);

    //split coffecients into evens and odds
    let evens: Vec<_> = coffecients.iter().step_by(2).cloned().collect();
    let odds: Vec<_> = coffecients.iter().skip(1).step_by(2).cloned().collect();
    let f0 = fft(&evens, &half_domain); //extract the even part,and do fft recursivly
    let f1 = fft(&odds, &half_domain); //extract the odd part,and do fft recursivly

    let mut left: Vec<FieldElement> = Vec::with_capacity(half_length);
    let mut right: Vec<FieldElement> = Vec::with_capacity(half_length);
    let factor = utils::get_factor_of_vec(domain);
    for i in 0..factor.len() {
        let l = f0[i] + factor[i] * f1[i];
        left.push(l);
        let r = f0[i] - factor[i] * f1[i];
        right.push(r);
    }
    right.reverse(); //in circle fft, we need let the right part do reverse.

    let mut o = vec![FieldElement::zero(); domain.len()];
    for i in 0..half_domain.len() {
        o[i] = left[i];
        o[i + half_domain.len()] = right[i];
    }
    return o;
}

// region: check poly by compute with basis
// compute the basis of Point(x,y)，that is bj(n),it's:1,y,x,yx,2x^2-1,y(2x^-1)......
fn fft_basis(point: Point, logn: u64) -> FieldElement {
    // let N=1<<logn-1;
    let mut jlist = utils::to_binary(logn);
    jlist.reverse(); //attention: jlist is the reverse of binary,such 4=001
    let y = point.y;
    let mut b = FieldElement::one();
    for i in 1..jlist.len() {
        if jlist[i] == 1 {
            b = b * point.v_n(i);
        }
    }
    if jlist.clone()[0] == 1 {
        b = b * y;
    }
    b
}

//evaluate at sigle point by coffecients (with circle_fft basis)
fn evaluate_by_fft_basis(coffecients: Vec<FieldElement>, point: Point) -> FieldElement {
    let mut out = FieldElement::zero();
    for i in 0..coffecients.len() {
        let j = i as u64;
        out = out + coffecients[i] * fft_basis(point, j)
    }
    out
}

fn check_poly(coffecients: Vec<FieldElement>, spcoset: Vec<Point>) -> Vec<FieldElement> {
    let mut results = Vec::new();
    for x in spcoset {
        results.push(evaluate_by_fft_basis(coffecients.clone(), x.clone())); // 这里是示例函数 f(x) = x * x
    }
    results
}
// endregion: MySection

#[test]
fn test() {
    // Generate standard position twin-cosets;
    let q = Point::new(FieldElement::new(24), FieldElement::new(18));
    let spcoset = Point::spcoset_by_Q(q, 3);
    for (index, value) in spcoset.iter().enumerate() {
        println!("spcoset:Index:{:?}, Value: {:?}", index, value);
    }
    println!("Generate Domain!");
    println!("------------------------------------");

    //example of ifft: from evaluations to coefficients
    let evaluations: Vec<FieldElement> = vec![
        FieldElement::new(1),
        FieldElement::new(2),
        FieldElement::new(4),
        FieldElement::new(8),
        FieldElement::new(16),
        FieldElement::new(1),
        FieldElement::new(2),
        FieldElement::new(4),
    ];
    // get coffecients by circle_ifft algorithm
    let coffecients = interpolate(&evaluations, &spcoset);
    let results = check_poly(coffecients.clone(), spcoset.clone()); //compute the evaluations with the binary polynomial
    assert_eq!(evaluations, results, "Values are not equal"); // check the result of circle_ifft algorithm is running correctly. If equal, then it's right.
    println!("get coffecients by circle_ifft successfully");
    println!("------------------------------------");

    // get evaluations by circle_fft algorithm
    let vals = evaluate(&coffecients, &spcoset);
    assert_eq!(evaluations, vals, "Values are not equal"); // check the result of circle_fft algorithm is running correctly. If equal, then it's right.
    println!("get evaluations by circle_fft successfully");
    println!("------------------------------------");

    //print the evaluations of fft
    println!("Now here are the evaluations!");
    for (index, evaluations) in vals.clone().iter().enumerate() {
        println!("vals:Index:{:?}, Value: {:?}", index, evaluations);
    }
    println!("------------------------------------");

    //print the coefficients
    println!("Now here are the coffecients!");
    for (index, value) in coffecients.clone().iter().enumerate() {
        println!("coffecients:Index:{:?}, Value: {:?}", index, value);
    }
    println!("------------------------------------");
}
