mod lib;

#[cfg(not(target_arch = "wasm32"))]
use std::fs;
#[cfg(not(target_arch = "wasm32"))]
use std::io::BufWriter;
#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;
#[cfg(not(target_arch = "wasm32"))]
use std::fs::File;
use text_io::scan;
use regex::Regex;
use lib::*;
use std::sync::{Arc, Mutex};

fn main() {
  let filename = "vertex_uvs";

  let contents = fs::read_to_string(filename).expect("Failed to read file");
  let lines: Vec<&str> = contents.split("\n").collect();

  let vertex_count: usize = lines[0].trim().parse().unwrap();
  let polygon_count: usize = lines[vertex_count+1].trim().parse().unwrap();
  let spline_count: usize = lines[vertex_count+polygon_count+2].trim().parse().unwrap();

  println!("{} vertices, {} polygons, {} splines", vertex_count, polygon_count, spline_count);

  let mut vertices: Vec<Vector> = Vec::new();
  let mut polygons: Vec<Vec<(usize, Vector2)>> = Vec::new();
  let mut splines: Vec<Vec<(Vector, Vector, Vector)>> = Vec::new();
  let mut spline_lengths: Vec<Vec<f64>> = Vec::new();

  for i in 0..vertex_count {
      let (x, y, z): (f64, f64, f64);
      scan!(lines[i+1].bytes() => "{} {} {}", x, y, z);
      vertices.push(Vector::new(x, y, z));
  }

  let loop_regex = Regex::new(r"(\d+) \((-?[\d]\.?[\d]+),\s*(-?[\d]\.?[\d]+)\)").unwrap();
  for i in 0..polygon_count {
      let mut poly = Vec::new();
      for capture in loop_regex.captures_iter(lines[i+2+vertex_count]) {
          let vertex_idx = capture[1].parse::<usize>().unwrap();
          let (uv_x, uv_y) = (capture[2].parse::<f64>().unwrap(), capture[3].parse::<f64>().unwrap());
          poly.push((vertex_idx, Vector2::new(uv_x, uv_y)));
      }
      polygons.push(poly);
  }

  let spline_regex = Regex::new(r"\(\((-?[\d]+\.?[\d]+),\s*(-?[\d]+\.?[\d]+),\s*(-?[\d]+\.?[\d]+)\),\s+\((-?[\d]+\.?[\d]+),\s*(-?[\d]+\.?[\d]+),\s*(-?[\d]+\.?[\d]+)\),\s+\((-?[\d]+\.?[\d]+),\s*(-?[\d]+\.?[\d]+),\s*(-?[\d]+\.?[\d]+)\)\)").unwrap();
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
              spline_length.push(spline_length[j-1] + bezier_test::bezier_distance(0., 1., &pts));
          }
      }
      spline_lengths.push(spline_length);
  }

  println!("{:?}", spline_lengths);

  // println!("{:?}", vertices);
  // println!("{:?}", polygons);
  // println!("{:?}", splines);

  // closest_quadratic_bezier_t_2d((-0.088, 0.845), vec![(0.226, 1.03), (0.5, 0.392), (0.34, 1.112)]);

  let (im_width, im_height) = (256, 256);

  let main_data = Arc::new(Mutex::new(vec![vec![vec![0_u16; 3]; im_width as usize]; im_height as usize]));

  let aux_data = Arc::new(Mutex::new(vec![vec![vec![0_f64; 3]; im_width as usize]; im_height as usize]));

  let mut test_tri = Vec::new();
  test_tri.push((0.3, 0.3));
  test_tri.push((0.3, 0.7));
  test_tri.push((0.6, 0.4));

  const STEP: usize = 1;
  let output = DataOutputType::LENGTH;

  let total_values = im_height * im_width;
  let completed_values = Arc::new(Mutex::new(0));
  let n_workers = 10;

  let pool = threadpool::ThreadPool::new(n_workers);

  let vertices = Arc::new(vertices);
  let polygons = Arc::new(polygons);
  let splines = Arc::new(splines);
  let spline_lengths = Arc::new(spline_lengths);
  let x_y_z_max = Arc::new(Mutex::new(0.0));
  
  for i in (0..im_height).step_by(STEP) {
      for j in (0..im_width).step_by(STEP) {
          // if tri_contains_point(&test_tri, (j as f32 / im_width as f32, i as f32 / im_height as f32)) {
          //     data[i as usize][j as usize] = [255, 255, 255].to_vec();
          // }
          let (vertices, polygons, splines, spline_lengths) = (vertices.clone(), polygons.clone(), splines.clone(), spline_lengths.clone());
          let x_y_z_max = x_y_z_max.clone();
          let (data, aux_data) = (main_data.clone(), aux_data.clone());
          let completed_values = completed_values.clone();
          pool.execute(move || {
              calculate_value_at_idx(i, j, &vertices, &polygons, &splines, &spline_lengths, x_y_z_max, im_width, im_height, data, aux_data, completed_values, STEP, output, total_values);
          });
      }
  }

  pool.join();

  extend_data(main_data.clone());
  extend_data(main_data.clone());

  #[cfg(not(target_arch = "wasm32"))]
  {
      let main_path = Path::new(r"./main_image.png");
      let main_file = File::create(main_path).unwrap();
      let ref mut main_w = BufWriter::new(main_file);

      let mut main_encoder = png::Encoder::new(main_w, im_width, im_height);
      main_encoder.set_color(png::ColorType::RGB);
      main_encoder.set_depth(png::BitDepth::Sixteen);
      let mut main_writer = main_encoder.write_header().unwrap();

      main_writer.write_image_data(&u16_vec_to_u8_vec(&flatten(flatten(Arc::try_unwrap(main_data).unwrap().into_inner().unwrap())))).unwrap();
      if output == DataOutputType::LENGTH {
          let aux_path = Path::new(r"./aux_image.png");
          let aux_file = File::create(aux_path).unwrap();
          let ref mut aux_w = BufWriter::new(aux_file);

          extend_data_f64(aux_data.clone());
          extend_data_f64(aux_data.clone());
          let mut aux_data = Arc::try_unwrap(aux_data).unwrap().into_inner().unwrap();
          for row in aux_data.iter_mut() {
              for px in row.iter_mut() {
                  for val in px.iter_mut() {
                      *val = 32768.0 * ((*val / *x_y_z_max.lock().unwrap()) + 1.0);
                  }
              }
          }
          
          let mut aux_encoder = png::Encoder::new(aux_w, im_width, im_height);
          aux_encoder.set_color(png::ColorType::RGB);
          aux_encoder.set_depth(png::BitDepth::Sixteen);
          let mut aux_writer = aux_encoder.write_header().unwrap();

          aux_writer.write_image_data(&u16_vec_to_u8_vec(&flatten(flatten(aux_data).into_iter().map(|x| vec![x[0] as u16, x[1] as u16, x[2] as u16]).collect()))).unwrap();
      }
  }
}
