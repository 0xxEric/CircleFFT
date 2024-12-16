use crate::field::FieldElement;
use crate::point::Point;

pub fn to_binary(n: u64) -> Vec<u64> {
    let binary_string = format!("{:b}", n);
    binary_string
        .chars()
        .map(|c| c.to_digit(10).unwrap() as u64)
        .collect()
}

pub fn get_halfdomain_of_point(domain: &Vec<Point>) -> Vec<FieldElement> {
    let half_length = domain.len() / 2;
    let half_domain: Vec<FieldElement> = domain
        .get(..half_length)
        .unwrap_or_default()
        .iter()
        .map(|p| p.x)
        .collect();
    return half_domain;
}

pub fn get_halfdomain_of_vec(domain: &Vec<FieldElement>) -> Vec<FieldElement> {
    let half_length = &domain.len() / 2;
    let half_domain: Vec<FieldElement> = domain
        .get(..half_length)
        .unwrap_or_default()
        .iter()
        .map(|p| (p.square().double() - FieldElement::one()))
        .collect();
    return half_domain;
}

pub fn get_factor_of_point(domain: &Vec<Point>) -> Vec<FieldElement> {
    let half_length = domain.len() / 2;
    let factor: Vec<FieldElement> = domain
        .get(..half_length)
        .unwrap_or_default()
        .iter()
        .map(|p| p.y)
        .collect();
    return factor;
}

pub fn get_factor_of_vec(domain: &Vec<FieldElement>) -> Vec<FieldElement> {
    let half_length = domain.len() / 2;
    let factor: Vec<FieldElement> = domain
        .get(..half_length)
        .unwrap_or_default()
        .iter()
        .map(|p| p.clone())
        .collect();
    return factor;
}

#[test]
fn test() {
    let a = to_binary(5);
    println!("a: {:?}", a);
    let b = to_binary(7);
    println!("b: {:?}", b);
}
