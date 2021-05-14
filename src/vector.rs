use std::ops;
use std::clone;
use std::cmp;

#[derive(Debug, Clone)]
pub struct Vector<T: ops::Add + clone::Clone> {
  values: Vec<T>
}

impl<T: ops::Add + clone::Clone> Vector<T> {
  pub fn new(values: &[T]) -> Vector<T> {
    Vector {
      values: values.to_vec(),
    }
  }
}

impl<T: ops::Add<Output = T> + clone::Clone> ops::Add for Vector<T> {
  type Output = Self;

  fn add(mut self, other: Self) -> Self::Output {
    let min = cmp::min(self.values.len(), other.values.len());
    for idx in 0..min {
      self.values[idx] = self.values[idx].clone() + other.values[idx].clone();
    }
    self
  }
}

impl<T: ops::Sub<Output = T> + clone::Clone + ops::Add> ops::Sub for Vector<T> {
  type Output = Self;

  fn sub(mut self, other: Self) -> Self::Output {
    let min = cmp::min(self.values.len(), other.values.len());
    for idx in 0..min {
      self.values[idx] = self.values[idx].clone() - other.values[idx].clone();
    }
    self
  }
}

impl<T: ops::Mul<Output = T> + clone::Clone + ops::Add> ops::Mul<T> for Vector<T> {
  type Output = Self;

  fn mul(mut self, other: T) -> Self::Output {
    for idx in 0..self.values.len() {
      self.values[idx] = self.values[idx].clone() * other.clone();
    }
    self
  }
}

impl<T: ops::Add + clone::Clone> ops::Index<usize> for Vector<T> {
  type Output = Option<T>;

  fn index(&self, idx: usize) -> &Option<T> {
    if self.values.len() > idx {
      return &Some(self.values[idx]);
    } else {
      return &None;
    }
  }
}