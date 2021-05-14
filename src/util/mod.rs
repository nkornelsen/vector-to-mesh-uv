use bezier_utils::BezierNum;

  type Vector2 = nalgebra::Vector2<BezierNum>;

  struct Line {
    slope: BezierNum,
    y_int: BezierNum,
    x: BezierNum,
    vertical: bool
  }

  impl Line {
    pub fn new(p1: Vector2, p2: Vector2) -> Line {
      let mut vertical = false;
      if p1[0] == p2[0] {
        vertical = true;
      }
      let mut slope = 0.;
      let mut x = 0.;
      let mut y_int = 0.;
      if !vertical {
        slope = (p2[1] - p1[1]) / (p2[0] - p1[0]);
        y_int = p1[1] - (slope * p1[0])
      } else {
        x = p1[0];
      }
      

      Line {
        slope,
        x,
        y_int,
        vertical
      }
    }

    pub fn above(&self, point: Vector2) -> bool {
      if !self.vertical {
        return point[1] > ((self.slope * point[0]) + self.y_int);
      } else {
        return point[0] > self.x;
      }
    }

    pub fn below(&self, point: Vector2) -> bool {
      if !self.vertical {
        return point[1] < ((self.slope * point[0]) + self.y_int);
      } else {
        return point[0] < self.x;
      }
    }
  }

  pub struct HalfPlane {
    line: Line,
    side: bool,
  }

  impl HalfPlane {
    pub fn new(p1: Vector2, p2: Vector2, side_point: Vector2) -> HalfPlane {
      let line = Line::new(p1, p2);
      let side = line.above(side_point);

      HalfPlane {
        line,
        side
      }
    }

    pub fn contains(&self, point: Vector2) -> bool {
      self.line.above(point) == self.side || self.line.below(point) != self.side
    }
  }