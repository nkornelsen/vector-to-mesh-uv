pub mod util;

use crate::util::HalfPlane;
use std::f32;
use std::io::{ Cursor, Seek, Read };
use std::sync::{Arc, Mutex};
use nalgebra::Vector3;
#[cfg(target_arch = "wasm32")]
use web_sys::console;

use std::io::BufWriter;
#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;
#[cfg(not(target_arch = "wasm32"))]
use std::fs::File;

#[cfg(not(target_arch = "wasm32"))]
use std::fs;
use text_io::scan;
use regex::Regex;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[derive(PartialEq, Copy, Clone)]
pub enum DataOutputType {
    BEZIER,
    DISTANCE,
    COORDS,
    LENGTH
}

pub type Vector = Vector3<f64>;
pub type Vector2 = nalgebra::Vector2<f64>;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn greet() {
    console::log_1(&"Hello from web-sys".into());
}

fn print(s: String) {
    #[cfg(target_arch = "wasm32")]
    console::log_1(&s.into());
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn generate_images(mesh_data: String, output_type: String) -> JsValue {
    let output = match output_type {
        "bezier" => DataOutputType::BEZIER,
        "distance" => DataOutputType::DISTANCE,
        "coords" => DataOutputType::COORDS,
        "length" => DataOutputType::LENGTH,
        _ => DataOutputType::BEZIER,
    };
    let (main, aux) = generate_images_internal(mesh_data, output);
    use js_sys::Uint8Array;
    let obj = js_sys::Object::new();
    js_sys::Reflect::set(&obj, &"main".into(), &Uint8Array::from(&main as &[u8])).unwrap();
    if let Some(a) = aux {
        js_sys::Reflect::set(&obj, &"aux".into(), &Uint8Array::from(&a as &[u8])).unwrap();
    } else {
        js_sys::Reflect::set(&obj, &"aux".into(), &JsValue::null()).unwrap();
    }
    JsValue::from(obj)
}

fn generate_images_internal(mesh_data: String, output_type: DataOutputType) -> (Vec<u8>, Option<Vec<u8>>) {
    // let filename = "vertex_uvs";
  
    // let contents = fs::read_to_string(filename).expect("Failed to read file");
    let lines: Vec<&str> = mesh_data.split("\n").collect();
  
    let vertex_count: usize = lines[0].trim().parse().unwrap();
    let polygon_count: usize = lines[vertex_count+1].trim().parse().unwrap();
    let spline_count: usize = lines[vertex_count+polygon_count+2].trim().parse().unwrap();
  
    print(format!("{} vertices, {} polygons, {} splines", vertex_count, polygon_count, spline_count));
  
    let mut vertices: Vec<Vector> = Vec::new();
    let mut polygons: Vec<Vec<(usize, Vector2)>> = Vec::new();
    let mut splines: Vec<Vec<(Vector, Vector, Vector)>> = Vec::new();
    let mut spline_lengths: Vec<Vec<f64>> = Vec::new();
  
    for i in 0..vertex_count {
        let (x, y, z): (f64, f64, f64);
        scan!(lines[i+1].bytes() => "{} {} {}", x, y, z);
        vertices.push(Vector::new(x, y, z));
    }
  
    let loop_regex = Regex::new(r"(\d+) \((-?[\d]*\.?[\d]*),\s*(-?[\d]*\.?[\d]*)\)").unwrap();
    for i in 0..polygon_count {
        let mut poly = Vec::new();
        for capture in loop_regex.captures_iter(lines[i+2+vertex_count]) {
            let vertex_idx = capture[1].parse::<usize>().unwrap();
            let (uv_x, uv_y) = (capture[2].parse::<f64>().unwrap(), capture[3].parse::<f64>().unwrap());
            poly.push((vertex_idx, Vector2::new(uv_x, uv_y)));
        }
        polygons.push(poly);
    }
  
    let spline_regex = Regex::new(r"\(\((-?[\d]*\.?[\d]*),\s*(-?[\d]*\.?[\d]*),\s*(-?[\d]*\.?[\d]*)\),\s+\((-?[\d]*\.?[\d]*),\s*(-?[\d]*\.?[\d]*),\s*(-?[\d]*\.?[\d]*)\),\s+\((-?[\d]*\.?[\d]*),\s*(-?[\d]*\.?[\d]*),\s*(-?[\d]*\.?[\d]*)\)\)").unwrap();
    for i in 0..spline_count {
        let mut spline = Vec::new();
        for capture in spline_regex.captures_iter(lines[i+3+vertex_count+polygon_count]) {
            let (p1, p2, p3) = 
                (Vector::new(capture[1].parse::<f64>().unwrap(), capture[2].parse::<f64>().unwrap(), capture[3].parse::<f64>().unwrap()), 
                Vector::new(capture[4].parse::<f64>().unwrap(), capture[5].parse::<f64>().unwrap(), capture[6].parse::<f64>().unwrap()), 
                Vector::new(capture[7].parse::<f64>().unwrap(), capture[8].parse::<f64>().unwrap(), capture[9].parse::<f64>().unwrap()));
            spline.push((p1, p2, p3));
        }
        splines.push(spline);
    }
  
    for i in 0..splines.len() {
        let mut spline_length = Vec::new();
        let spline = &splines[i];
        spline_length.push(0.0);
        if spline.len() > 1 {
            for j in 1..spline.len() {
                let pts = vec![spline[j-1].1, spline[j-1].2, spline[j].0, spline[j].1];
                spline_length.push(spline_length[j-1] + bezier_utils::bezier_distance(0., 1., &pts));
            }
        }
        spline_lengths.push(spline_length);
    }
  
    print(format!("{:?}", spline_lengths));
  
    print(format!("{:?}", vertices));
    print(format!("{:?}", polygons));
    print(format!("{:?}", splines));
  
    // closest_quadratic_bezier_t_2d((-0.088, 0.845), vec![(0.226, 1.03), (0.5, 0.392), (0.34, 1.112)]);
  
    let (im_width, im_height) = (256, 256);
  
    let mut main_data = vec![vec![vec![0_u16; 3]; im_width as usize]; im_height as usize];
  
    let mut aux_data = vec![vec![vec![0_f64; 3]; im_width as usize]; im_height as usize];

    let mut write_mask = vec![vec![false; im_width as usize]; im_height as usize];
  
    let mut test_tri = Vec::new();
    test_tri.push((0.3, 0.3));
    test_tri.push((0.3, 0.7));
    test_tri.push((0.6, 0.4));
  
    const STEP: usize = 1;
    let output = output_type;
  
    let total_values = im_height * im_width;
    let mut completed_values = 0;
  
    print(format!("Ready for processing"));
  
    let vertices = vertices;
    let polygons = polygons;
    let splines = splines;
    let spline_lengths = spline_lengths;
    let mut x_y_z_max = 0.0;

    
    for i in (0..im_height).step_by(STEP) {
        for j in (0..im_width).step_by(STEP) {
            // if tri_contains_point(&test_tri, (j as f32 / im_width as f32, i as f32 / im_height as f32)) {
            //     data[i as usize][j as usize] = [255, 255, 255].to_vec();
            // }
            calculate_value_at_idx(i, j, &vertices, &polygons, &splines, &spline_lengths, &mut x_y_z_max, im_width, im_height, &mut main_data, &mut aux_data, &mut write_mask, &mut completed_values, STEP, output, total_values);
        }
    }
  
    // pool.join();
  
    extend_data(&mut main_data, &mut write_mask);
    extend_data(&mut main_data, &mut write_mask);
    {
        let mut main_image_data = Vec::new();

        {
            let mut main_encoder = png::Encoder::new(&mut main_image_data, im_width, im_height);
            main_encoder.set_color(png::ColorType::RGB);
            main_encoder.set_depth(png::BitDepth::Sixteen);
            let mut main_writer = main_encoder.write_header().unwrap();
            
            main_writer.write_image_data(&u16_vec_to_u8_vec(&flatten(flatten(main_data)))).unwrap();
        }
        

        let aux_output = if output == DataOutputType::LENGTH {
            let mut aux_image_data = Vec::new();
  
            extend_data_f64(&mut aux_data, &mut write_mask);
            extend_data_f64(&mut aux_data, &mut write_mask);
            for row in aux_data.iter_mut() {
                for px in row.iter_mut() {
                    for val in px.iter_mut() {
                        *val = 32768.0 * ((*val / x_y_z_max) + 1.0);
                    }
                }
            }
            
            {
                let mut aux_encoder = png::Encoder::new(&mut aux_image_data, im_width, im_height);
                aux_encoder.set_color(png::ColorType::RGB);
                aux_encoder.set_depth(png::BitDepth::Sixteen);
                let mut aux_writer = aux_encoder.write_header().unwrap();
    
                aux_writer.write_image_data(&u16_vec_to_u8_vec(&flatten(flatten(aux_data).into_iter().map(|x| vec![x[0] as u16, x[1] as u16, x[2] as u16]).collect()))).unwrap();
            }

            Some(aux_image_data)
        } else {
            None
        };
        (main_image_data, aux_output)
    }
}

pub fn u16_vec_to_u8_vec(data: &Vec<u16>) -> Vec<u8> {
    let mut output = vec![0_u8; data.len() * 2];
    for i in 0..data.len() {
        output[2*i] = (data[i] >> 8) as u8;
        output[2*i+1] = (data[i] & 0x00FF) as u8;
    }
    output
}

pub fn extend_data(data: &mut Vec<Vec<Vec<u16>>>, write_mask: &mut Vec<Vec<bool>>) {
    let mut overlay_data = vec![vec![vec![0_u16; 3]; data[0].len()]; data.len()];

    for i in 0..data.len() as isize {
        for j in 0..data.len() as isize {
            let (ui, uj) = (i as usize, j as usize);
            // if data[ui][uj][0] != 0 || data[ui][uj][1] != 0 || data[ui][uj][2] != 0 {
            if write_mask[ui][uj] {
                let offsets = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)];
                for o in offsets.iter() {
                    let (oi, oj) = (i+o.0, j+o.1);
                    if oi < data.len() as isize && oj < data[0].len() as isize && oi >= 0 && oj > 0 {
                        let (oi, oj) = (oi as usize, oj as usize);
                        if !write_mask[oi][oj] {
                            overlay_data[oi][oj][0] = data[ui][uj][0];
                            overlay_data[oi][oj][1] = data[ui][uj][1];
                            overlay_data[oi][oj][2] = data[ui][uj][2];
                            write_mask[oi][oj] = true;
                        }
                    }
                }
            }
        }
    }

    for i in 0..data.len() {
        for j in 0..data.len() {
            for k in 0..3 {
                data[i][j][k] += overlay_data[i][j][k];
            }
        }
    }
}

pub fn extend_data_f64(data: &mut Vec<Vec<Vec<f64>>>, write_mask: &Vec<Vec<bool>>) {
    let mut overlay_data = vec![vec![vec![0_f64; 3]; data[0].len()]; data.len()];

    for i in 0..data.len() as isize {
        for j in 0..data.len() as isize {
            let (ui, uj) = (i as usize, j as usize);
            if write_mask[ui][uj] {
                let offsets = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)];
                for o in offsets.iter() {
                    let (oi, oj) = (i+o.0, j+o.1);
                    if oi < data.len() as isize && oj < data[0].len() as isize && oi >= 0 && oj > 0 {
                        let (oi, oj) = (oi as usize, oj as usize);
                        if !write_mask[oi][oj] {
                            overlay_data[oi][oj][0] = data[ui][uj][0];
                            overlay_data[oi][oj][1] = data[ui][uj][1];
                            overlay_data[oi][oj][2] = data[ui][uj][2];
                        }
                    }
                }
            }
        }
    }

    for i in 0..data.len() {
        for j in 0..data.len() {
            for k in 0..3 {
                data[i][j][k] += overlay_data[i][j][k];
            }
        }
    }
}

pub fn calculate_value_at_idx(i: u32, j: u32, 
        vertices: &Vec<Vector>, 
        polygons: &Vec<Vec<(usize, Vector2)>>, 
        splines: &Vec<Vec<(Vector, Vector, Vector)>>,
        spline_lengths: &Vec<Vec<f64>>,
        x_y_z_max: &mut f64,
        im_width: u32, im_height: u32, 
        data: &mut Vec<Vec<Vec<u16>>>,
        aux_data: &mut Vec<Vec<Vec<f64>>>,
        write_mask: &mut Vec<Vec<bool>>,
        _completed_values: &mut isize, 
        step: usize, 
        output: DataOutputType, 
        total_values: u32) {
    for k in 0..polygons.len() {
        let point = Vector2::new((j as f64 + 0.5) / im_width as f64, (i as f64 + 0.5) / im_height as f64);
        if tri_contains_point(&vec![polygons[k][0].1, polygons[k][1].1, polygons[k][2].1], point) {
            // data[i as usize][j as usize] = [k as u8 * 30, k as u8 * 30, k as u8 * 30].to_vec();
            let barycentric = cartesian_to_barycentric(&polygons[k], point);
            let mut triangle = Vec::new();
            for l in 0..polygons[k].len() {
                triangle.push(vertices[polygons[k][l].0]);
            }
            // println!("Using triangle: {:?}", triangle);
            let cartesian = barycentric_to_cartesian(&triangle, barycentric);

            let best_spline = closest_spline(&splines, cartesian);
            if step > 1 {
                println!("Triangle: {:?}", polygons[k]);
                println!("{:?} -> {:?} -> {:?}", point, barycentric, cartesian);
                println!("Using spline: {:?}", best_spline);
            }
            let mut derivative = bezier_derivative(&splines, best_spline);
            // let mut derivative = bezier_utils::point_derivative(best_spline.0, );
            for d in derivative.iter_mut() {
                *d = (*d + 1.0) * 32768_f64;
            }
            // println!("Derivative: {:?}", derivative);
            // data[ im_height as usize - 1 - i as usize][j as usize] = [((derivative[0] + 1.0) * 128_f32) as u8, ((derivative[1] + 1.0) * 128_f32) as u8, ((derivative[2] + 1.0) * 128_f32) as u8 ].to_vec();
            {
                write_mask[im_height as usize - 1 - i as usize][j as usize] = true;
                if output == DataOutputType::DISTANCE {
                    data[im_height as usize - 1 - i as usize][j as usize] = [(best_spline.3 * 32768_f64) as u16, (best_spline.3 * 32768_f64) as u16, (best_spline.3 * 32768_f64) as u16].to_vec();
                } else if output == DataOutputType::BEZIER {
                    data[im_height as usize - 1 - i as usize][j as usize] = [derivative[0] as u16, derivative[1] as u16, derivative[2] as u16].to_vec();
                } else if output == DataOutputType::COORDS {
                    let fac = 32768_f64;
                    data[im_height as usize - 1 - i as usize][j as usize] = [((cartesian[0] + 1_f64) * fac) as u16, ((cartesian[1] + 1_f64) * fac) as u16, ((cartesian[2] + 1_f64) * fac) as u16].to_vec();
                } else if output == DataOutputType::LENGTH {
                    let (a, b, t, _) = best_spline;
                    let pts = vec![splines[a][b].1, splines[a][b].2, splines[a][b+1].0, splines[a][b+1].1];
                    let length_on_current_curve = 65535.0 * (spline_lengths[a][b] + bezier_utils::bezier_distance(0., t, &pts)) / spline_lengths[a][spline_lengths[a].len() - 1];
                    
                    // println!("length: {:?}", length_on_current_curve);
                    data[im_height as usize - 1 - i as usize][j as usize] = [length_on_current_curve as u16, length_on_current_curve as u16, length_on_current_curve as u16].to_vec();
                    // also store offset
                    let offset = cartesian - bezier_utils::point_location(best_spline.2, &pts);
                    (*aux_data)[im_height as usize - 1 - i as usize][j as usize] = [offset[0], offset[1], offset[2]].to_vec();
                    for i in 0..3 {
                        if offset[i].abs() > *x_y_z_max {
                            *x_y_z_max = offset[i].abs();
                        }
                    }
                }
                // data[ im_height as usize - 1 - i as usize][j as usize] = [128_u8, 128_u8, 128_u8].to_vec();
                break;
            }
        }
    }
    {
        *_completed_values += 1;

        if *_completed_values % (((total_values + 100) / 100) as isize) == 0 {
            print(format!("{:.2}% done", 100_f32 * (*_completed_values as f32) / total_values as f32));
        }
    }
}

fn cartesian_to_barycentric(triangle: &Vec<(usize, Vector2)>, point: Vector2) -> Vector {
    let (x1, y1) = (triangle[0].1[0], triangle[0].1[1]);
    let (x2, y2) = (triangle[1].1[0], triangle[1].1[1]);
    let (x3, y3) = (triangle[2].1[0], triangle[2].1[1]);
    let (x, y) = (point[0], point[1]);

    let det = ((y2 - y3) * (x1 - x3)) + ((x3 - x2) * (y1 - y3));
    let l1 = (((y2 - y3) * (x - x3)) + ((x3 - x2) * (y - y3))) / det;
    let l2 = (((y3 - y1) * (x - x3)) + ((x1 - x3) * (y - y3))) / det;
    let l3 = 1_f64 - l1 - l2;

    Vector::new(l1, l2, l3)
}

fn barycentric_to_cartesian(triangle: &Vec<Vector>, point: Vector) -> Vector {
    let (mut l1, mut l2, mut l3) = (point[0], point[1], point[2]);
    if l1 + l2 + l3 != 1_f64 {
        // println!("Normalizing");
        let s = l1 + l2 + l3;
        l1 /= s;
        l2 /= s;
        l3 /= s;
    }
    let (x1, y1, z1) = (triangle[0][0], triangle[0][1], triangle[0][2]);
    let (x2, y2, z2) = (triangle[1][0], triangle[1][1], triangle[1][2]);
    let (x3, y3, z3) = (triangle[2][0], triangle[2][1], triangle[2][2]);

    // println!("Finding {:?} in {:?} {:?} {:?}", (l1, l2, l3), (x1, y1, z1), (x2, y2, z2), (x3, y3, z3));

    let x = (l1*x1) + (l2*x2) + (l3*x3);
    let y = (l1*y1) + (l2*y2) + (l3*y3);
    let z = (l1*z1) + (l2*z2) + (l3*z3);

    Vector::new(x, y, z)
}

fn tri_contains_point(triangle: &Vec<Vector2>, point: Vector2) -> bool {
    let hp1 = HalfPlane::new(triangle[0], triangle[1], triangle[2]);
    let hp2 = HalfPlane::new(triangle[2], triangle[0], triangle[1]);
    let hp3 = HalfPlane::new(triangle[1], triangle[2], triangle[0]);

    hp1.contains(point) && hp2.contains(point) && hp3.contains(point)
}

pub fn flatten<T>(nested: Vec<Vec<T>>) -> Vec<T> {
    nested.into_iter().flatten().collect()
}

fn closest_spline(splines: &Vec<Vec<(Vector, Vector, Vector)>>, point: Vector) -> (usize, usize, f64, f64) {
    let mut best = (0, 0, 0_f64, -1_f64);
    for i in 0..splines.len() {
        let spline = &splines[i];
        for j in 0..spline.len()-1 {
            let p0 = spline[j].1;
            let p1 = spline[j].2;
            let p2 = spline[j+1].0;
            let p3 = spline[j+1].1;

            // let cur = closest_quadratic_bezier_t_3d(point, vec![p0, p1, p2]);
            let cur = bezier_utils::distance_to_curve(point, &vec![p0, p1, p2, p3]);
            if best.3 == -1_f64 || cur.1 < best.3 {
                best = (i, j, cur.0, cur.1);
            }
        }
    }
    best
}

fn bezier_derivative(splines: &Vec<Vec<(Vector, Vector, Vector)>>, param: (usize, usize, f64, f64)) -> Vector {
    // use nalgebra::Vector3;
    // let p0 = Vector3::new(splines[param.0][param.1].1.0, splines[param.0][param.1].1.1, splines[param.0][param.1].1.2);
    // let p1 = Vector3::new(splines[param.0][param.1].2.0, splines[param.0][param.1].2.1, splines[param.0][param.1].2.2);
    // let p2 = Vector3::new(splines[param.0][param.1+1].1.0, splines[param.0][param.1+1].1.1, splines[param.0][param.1+1].1.2);
    // let t = param.2;
    // // println!("Derivative calculating for {:?} {:?} {:?} @ {:?}", p0, p1, p2, t);

    // let mut deriv = 2_f32*(1_f32 - t)*(p1 - p0) + 2_f32*t*(p2 - p1);
    // deriv = deriv.normalize() + Vector3::new(1_f32, 1_f32, 1_f32);
    // vec![deriv.x, deriv.y, deriv.z]

    let p0 = splines[param.0][param.1].1;
    let p1 = splines[param.0][param.1].2;
    let p2 = splines[param.0][param.1+1].0;
    let p3 = splines[param.0][param.1+1].1;
    
    let t = param.2;
    let deriv = bezier_utils::point_derivative(t, &vec![p0, p1, p2, p3]).normalize();

    deriv
}