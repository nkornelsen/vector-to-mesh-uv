use std::collections::LinkedList;

type BezierNum = f64;
type Vector<T> = nalgebra::Vector3<T>;

const ZERO_TOLERANCE: BezierNum = 0.00000000001;
const PRINT_DEBUG: bool = false;

fn print_polynomials(a: &Vec<Vec<BezierNum>>, n: usize) {
    for i in 0..(n+1) {
        print!("S_{}(u) = ", i);
        for j in 0..(n-i + 1) {
            if n - i - j > 0 {
                print!("{:?}u^{}", a[j][i], n - i - j);
            } else {
                print!("{:?}", a[j][i]);
            }
            if (j as isize) < ((n - i) as isize) {
                print!(" + ");
            }
        }
        println!("");
    }
}

fn calculate_polynomial_value(t: BezierNum, idx: usize, coefficients: &Vec<Vec<BezierNum>>, degree: usize) -> BezierNum {
    let mut value = 0.;
    for j in 0..(degree - idx + 1) {
        let pow = degree - idx - j;
        value += coefficients[j][idx] * t.powi(pow as i32);
    }
    value
}

fn sturm_sequence_at(t: BezierNum, coefficients: &Vec<Vec<BezierNum>>, degree: usize) -> Vec<i32> {
    let mut output = vec![0; 6];

    for i in 0..(degree+1) {
        let value = calculate_polynomial_value(t, i, coefficients, degree);

        output[i] = {
            if value.abs() <= ZERO_TOLERANCE {
                0
            } else if value > 0. {
                1
            } else {
                -1
            }
        };
    }

    output
}

fn sign_changes(mut signs: Vec<i32>) -> i32 {
    signs = signs.into_iter().filter(|e| *e != 0).collect();

    let mut changes = 0;
    for i in 1..signs.len() {
        if signs[i-1] != signs[i] {
            changes += 1;
        }
    }
    
    changes
}

fn print_sign_sequence(signs: &Vec<i32>) {
    for i in 0..signs.len() {
        print!("{} ", if signs[i] > 0 {"+"} else if signs[i] < 0 {"-"} else {"0"});
    }
    println!("");
}

fn sturm_count_roots(a: BezierNum, b: BezierNum, coefficients: &Vec<Vec<BezierNum>>, degree: usize) -> i32 {
    sign_changes(sturm_sequence_at(a, coefficients, degree)) - sign_changes(sturm_sequence_at(b, coefficients, degree))
}

fn locate_root_intervals(a: BezierNum, b: BezierNum, coefficients: &Vec<Vec<BezierNum>>, degree: usize) -> Vec<(BezierNum, BezierNum)> {
    let mut intervals = vec![];

    let check_interval_signs = |i: (BezierNum, BezierNum)| {
        calculate_polynomial_value(i.0, 0, coefficients, degree) <= 0. &&
        calculate_polynomial_value(i.1, 0, coefficients, degree) > 0.
    };

    let initial_roots = sturm_count_roots(a, b, coefficients, degree);
    if initial_roots == 0 {
        intervals
    } else if initial_roots == 1 {
        if check_interval_signs((a, b)) {
            intervals.push((a, b));
        }
        intervals
    } else {
        let mut roots_queue = LinkedList::new();
        roots_queue.push_back((a, b));

        while !roots_queue.is_empty() {
            let interval = roots_queue.pop_front().unwrap();
            let midpoint = (interval.0 + interval.1) / 2.0;
            
            let i1 = (interval.0, midpoint);
            let i2 = (midpoint, interval.1);

            let mut push_if_valid = |i: (BezierNum, BezierNum)| {
                let roots = sturm_count_roots(i.0, i.1, coefficients, degree);
                if roots == 1 {
                    if calculate_polynomial_value(i.1, 0, coefficients, degree) == 0. {
                        intervals.push((i.1, i.1));
                    } else if check_interval_signs(i) {
                        intervals.push(i);
                    }
                } else if roots > 1 {
                    roots_queue.push_back(i);
                }
            };

            push_if_valid(i1);
            push_if_valid(i2);
        }
        intervals
    }
}

fn compute_root_on_interval(int: (BezierNum, BezierNum), coefficients: &Vec<f64>, degree: usize) -> f64 {
    let tolerance = 0.0000001;
    let int = (int.0, int.1);
    // println!("Coefficients: {:?}", coefficients);
    let mut a = int.0;
    let mut b = int.1;
    let mut r = (a + b) * 0.5;
    use std::convert::TryInto;

    while b - a > tolerance {
        let m = (a + b) * 0.5;
        let g_m = {
            let mut value = 0.;
            for i in 0..coefficients.len() {
                value += coefficients[i] * m.powi((degree - i).try_into().unwrap());
            }
            value
        };
        // println!("g({}) = {:?}", m, g_m);
        if g_m == 0. {
            r = m;
            break;
        }
        if g_m < 0. {
            a = m;
        } else {
            b = m;
        }
        r = (a + b) * 0.5;
    }
    r
}

fn roots_from_intervals(intervals: Vec<(BezierNum, BezierNum)>, coefficients: &Vec<f64>, degree: usize) -> Vec<f64> {
    intervals.iter().map(|i| {
        compute_root_on_interval(*i, coefficients, degree)
    }).collect()
}

fn choose_best_root(roots: &Vec<f64>, point: Vector<BezierNum>, points: &Vec<Vector<BezierNum>>) -> (f64, f64) {
    // let point = (point[0], point[1]);
    let dist_sq = |r: f64| {
        let p = point_location(r, points);
        (p - point).magnitude_squared()
    };

    let mut best = (-1., 0.);

    for &root in roots {
        if best.0 == -1. || dist_sq(root) < best.1 {
            best = (root, dist_sq(root));
        }
    }
    (best.0, best.1.sqrt())
}

pub fn point_location(t: f64, points: &Vec<Vector<BezierNum>>) -> Vector<BezierNum> {
    // let points: Vec<(f64, f64)> = points.iter().map(|x| (x[0], x[1])).collect();

    (1.0-t).powi(3)*points[0] + 3.0*(1.0-t).powi(2)*t*points[1] + 3.0*(1.0-t)*t.powi(2)*points[2] + t.powi(3)*points[3]
}

pub fn point_derivative(t: f64, points: &Vec<Vector<BezierNum>>) -> Vector<BezierNum> {
    3.0*(1.0-t).powi(2)*(points[1]-points[0]) + 6.0*(1.0-t)*t*(points[2]-points[1]) + 3.0*t.powi(2)*(points[3]-points[2])
}

fn print_ratios(ratios: &Vec<BezierNum>) {
    for r in ratios {
        print!("{} ", r);
    }
    println!("");
}

pub fn distance_to_curve(p: Vector<BezierNum>, points: &Vec<Vector<BezierNum>>) -> (BezierNum, BezierNum) {
    let n = -points[0] + points[1]*3. - points[2]*3. + points[3];
    let r = points[0]*3. - points[1]*6. + points[2]*3.;
    let s = points[0]*-3. + points[1]*3.;
    let v = points[0];

    let j = n*3.;
    let k = r*2.;

    // let p = Vector::new(BigRational::from_float(0.5).unwrap(), BigRational::from_float(0.25).unwrap());
    let b = vec![n.dot(&j), n.dot(&k)+r.dot(&j), n.dot(&s)+r.dot(&k)+s.dot(&j), -p.dot(&j)+r.dot(&s)+s.dot(&k)+v.dot(&j), -p.dot(&k)+s.dot(&s)+v.dot(&k), -p.dot(&s)+v.dot(&s)];
    // let b_float: Vec<f64> = b.iter().map(|x| x).collect();

    // let b_float = vec![1000., -1500., 850., -225., 27.4, -1.2];
    if PRINT_DEBUG {
        println!("{:?}", b);
    }
    let n = 5;

    let mut a: Vec<Vec<BezierNum>> = vec![vec![0.; n+1]; n+2];
    let mut t: Vec<BezierNum> = vec![0.; n+1];
    let mut m: Vec<BezierNum> = vec![0.; n+1];
    
    for j in 0..n+1 {
        a[j][0] = b[j];
        a[j][1] = a[j][0] * (n - j) as f64;
    }
    t[2] = a[0][0] / a[0][1];
    m[2] = (a[1][0] - t[2] * a[1][1]) / a[0][1];

    for i in 2..n+1 {
        for j in 0..n {
            a[j][i] = -(a[j+2][i-2] - m[i]*a[j+1][i-1] - t[i]*a[j+2][i-1]);
        }
        if i < n {
            let i = i+1;
            if a[0][i-1] != 0. {
                t[i] = a[0][i-2] / a[0][i-1];
                m[i] = (a[1][i-2] - (t[i]*a[1][i-1])) / a[0][i-1];
            }
        }
    }

    if PRINT_DEBUG {
        println!("a: {:?}\nt: {:?}\nm: {:?}", a, t, m);
        print_polynomials(&a, n);
        print_sign_sequence(&sturm_sequence_at(0., &a, n));
        print_sign_sequence(&sturm_sequence_at(1., &a, n));
    }
    let mut test_roots = roots_from_intervals(locate_root_intervals(0., 1., &a, n), &b, n);
    if PRINT_DEBUG {
        println!("Roots on (0, 1]: {:?}", test_roots);
    }

    test_roots.push(1.0);
    test_roots.push(0.0);

    choose_best_root(&test_roots, p, &points)
}

fn flatten<T>(nested: Vec<Vec<T>>) -> Vec<T> {
    nested.into_iter().flatten().collect()
}

const DISTANCE_TOLERANCE: BezierNum = 0.0000001;

fn bezier_distance_recurse(a: BezierNum, b: BezierNum, points: &Vec<Vector<BezierNum>>, direct_dist: BezierNum) -> BezierNum {
    let midpoint = (a + b) / 2.0;
    let s1_dist = line_distance_on_curve(a, midpoint, points);
    let s2_dist = line_distance_on_curve(midpoint, b, points);
    let divided_dist = s1_dist + s2_dist;

    if (direct_dist - divided_dist).abs() < DISTANCE_TOLERANCE {
        return direct_dist;
    } else {
        return bezier_distance_recurse(a, midpoint, points, s1_dist) + bezier_distance_recurse(midpoint, b, points, s2_dist);
    }
}

fn line_distance_on_curve(a: BezierNum, b: BezierNum, points: &Vec<Vector<BezierNum>>) -> BezierNum {
    let p1 = point_location(a, points);
    let p2 = point_location(b, points);
    // println!("{:?} -> {:?} = {}", p1, p2, (p1 - p2).magnitude());
    (p1 - p2).magnitude()
}

pub fn bezier_distance(a: BezierNum, b: BezierNum, points: &Vec<Vector<BezierNum>>) -> BezierNum {
    let direct_dist = line_distance_on_curve(a, b, points);
    let midpoint = (a + b) / 2.0;
    let s1_dist = line_distance_on_curve(a, midpoint, points);
    let s2_dist = line_distance_on_curve(midpoint, b, points);
    let divided_dist = s1_dist + s2_dist;

    if (direct_dist - divided_dist).abs() < DISTANCE_TOLERANCE {
        return direct_dist;
    } else {
        return bezier_distance_recurse(a, midpoint, points, s1_dist) + bezier_distance_recurse(midpoint, b, points, s2_dist);
    }
}