use vector_to_mesh_uv::*;
use std::fs::File;
use std::path::Path;
use std::io::{Read, Write};
use clap::{Arg, App, SubCommand};

fn main() {
    let matches = App::new("vector-to-mesh-uv")
        .arg(Arg::with_name("input")
            .short("i")
            .long("input")
            .value_name("INPUT")
            .help("The input mesh data")
            .takes_value(true))
        .arg(Arg::with_name("main_out")
            .short("o")
            .long("main_out")
            .value_name("MAIN_OUT")
            .help("Output file for main png")
            .takes_value(true))
        .arg(Arg::with_name("aux_out")
            .short("aux")
            .long("aux_out")
            .value_name("AUXILIARY_OUT")
            .help("Output file for secondary png")
            .takes_value(true))
        .arg(Arg::with_name("output_type")
            .short("t")
            .long("output_type")
            .value_name("OUTPUT_TYPE")
            .help("Type of output to generate")
            .takes_value(true))
        .arg(Arg::with_name("width")
            .short("w")
            .long("width")
            .value_name("WIDTH")
            .help("Width of output")
            .takes_value(true))
        .arg(Arg::with_name("height")
            .short("h")
            .long("height")
            .value_name("HEIGHT")
            .help("Height of output")
            .takes_value(true))
        .arg(Arg::with_name("multithreaded")
            .short("mt")
            .long("multithread")
            .value_name("MULTITHREADING")
            .help("Use multiple threads")
            .takes_value(true))
        .arg(Arg::with_name("workers")
            .short("n")
            .long("workers")
            .value_name("WORKERS")
            .help("Number of workers (if multithreading)")
            .takes_value(true)
    ).get_matches();

//     let (main, aux) = generate_images_internal_multithreaded(String::from("4
// -1 -1 1
// 1 -1 1
// -1 1 1
// 1 1 1
// 2
// 1 (0, 1) 2 (1, 0) 0 (1, 1) 
// 1 (0, 1) 3 (0, 0) 2 (1, 0) 
// 1
// ((-1.1589628458023071, -1.6762653589248657, 1), (-0.6589628458023071, -1.1762653589248657, 1), (2.4685609340667725, 1.3283957242965698, 1)) ((-3.308171272277832, 1.3792165517807007, 1), (1.7162351608276367, -1.2275062799453735, 1), (2.963348150253296, 0.1608094573020935, 1))
// "), DataOutputType::BEZIER, 512, 255, 10);

//     let main_path = Path::new("./main_image.png");
//     let mut main_file = File::create(main_path).unwrap();
//     main_file.write(&main).unwrap();

//     if let Some(aux) = aux {
//         let aux_path = Path::new("./aux_image.png");
//         let mut aux_file = File::create(aux_path).unwrap();
//         aux_file.write(&aux).unwrap();
//     }
    let mut input_file = File::open(matches.value_of("input").unwrap()).unwrap();
    let mut input_data = String::new();
    input_file.read_to_string(&mut input_data).unwrap();

    let width: u32 = matches.value_of("width").unwrap_or("255").parse().unwrap();
    let height: u32 = matches.value_of("height").unwrap_or("255").parse().unwrap();

    let multithreaded: bool = matches.value_of("multithreaded").unwrap_or("true").parse().unwrap();
    let workers: usize = matches.value_of("workers").unwrap_or("4").parse().unwrap();

    let output_type = matches.value_of("output_type").unwrap_or("bezier");
    let output = match output_type {
        "bezier" => DataOutputType::BEZIER,
        "distance" => DataOutputType::DISTANCE,
        "coords" => DataOutputType::COORDS,
        "length" => DataOutputType::LENGTH,
        _ => DataOutputType::BEZIER,
    };

    let (main, aux) = if multithreaded {
        generate_images_internal_multithreaded(input_data, output, width, height, workers)
    } else {
        generate_images_internal(input_data, output, width, height)
    };

    let main_path = Path::new(matches.value_of("main_out").unwrap_or("./main_image.png"));
    let mut main_file = File::create(main_path).unwrap();
    main_file.write(&main).unwrap();

    if let Some(aux) = aux {
        let aux_path = Path::new(matches.value_of("aux_out").unwrap_or("./aux_image.png"));
        let mut aux_file = File::create(aux_path).unwrap();
        aux_file.write(&aux).unwrap();
    }

}