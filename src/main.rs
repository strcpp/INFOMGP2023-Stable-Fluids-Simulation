use macroquad::prelude::*;
mod fluid;

fn window_conf() -> Conf {
    Conf {
        window_title: "Fluid thing?".to_owned(),
        window_width: 500,
        window_height: 500, 
        window_resizable: false,
        ..Default::default()
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    let WIDTH = screen_width() as usize;
    let HEIGHT = screen_height() as usize;

    let mut pixels = vec![0; WIDTH * HEIGHT * 4];

    let mut fluid_grid = fluid::FluidGrid::new(WIDTH as usize, HEIGHT as usize);
    let dt = 0.1;;
    let fluid = fluid::Fluid::new(WIDTH as usize, HEIGHT as usize, dt, 0.0, 0.0);

    let texture = Texture2D::from_rgba8(WIDTH as u16, HEIGHT as u16, &pixels);
    
    let mut image =  Image {
        bytes: pixels,
        width: WIDTH as u16,
        height: HEIGHT as u16   
    };

    
    let mut last_mouse_position: Option<(f32, f32)> = None;
    let mut pvx = 0.0;
    let mut pvy = 0.0;

    let mut current_color = Color::new(1.0, 0.0, 0.0, 1.0);

    let mut time_counter = 0.0;
    let radius = 10.0;
    let num_points = 30;
    let center_x = WIDTH as f32 / 2.0;
    let center_y = HEIGHT as f32 - radius;
    let mut vort = true;
    let mut viz_temp = false;

    loop {
        let fps = get_fps();

        let mut prev_x = 0;
        for i in 0..num_points {
            let angle = i as f32 * 2.0 * std::f32::consts::PI / num_points as f32;
            let random_velocity_x = rand::gen_range(-1.0, 1.0);
            let random_velocity_y = rand::gen_range(-1.0, 1.0);
            let x = prev_x + WIDTH / (num_points + 1);
            let y = HEIGHT - 1;
            let random_density = rand::gen_range(0.0, 7.0);
            prev_x = x;

            fluid_grid.add_velocity(x as usize, y as usize, random_velocity_x as f64, random_velocity_y as f64);
            fluid_grid.add_density(x as usize, y as usize, random_density, current_color);
        }
        
        clear_background(BLACK);
        fluid_grid.apply_gravity(-0.1 * (1.0/fps as f64));


        if is_key_pressed(KeyCode::Escape) {
            vort = !vort;
        }

        if(is_key_pressed(KeyCode::Q)) {
            viz_temp = !viz_temp;
        }

        if(is_key_pressed(KeyCode::Space)) {
            current_color = Color::new(
                rand::gen_range(0.0,1.0),
                rand::gen_range(0.0,1.0),
                rand::gen_range(0.0,1.0),
                1.0,
            );
        }

        fluid.step(&mut fluid_grid, vort);
        fluid_grid.render_densities(&mut image, &current_color, viz_temp);

        if is_mouse_button_down(MouseButton::Left) {
            let current_mouse_position = mouse_position();

            if let Some(last_position) = last_mouse_position {
                let vx = (current_mouse_position.0 - last_position.0) * 2.0;
                let vy = (current_mouse_position.1 - last_position.1) * 2.0;

                fluid_grid.add_density(current_mouse_position.0 as usize, current_mouse_position.1 as usize, 100.0, current_color);
                fluid_grid.add_velocity(current_mouse_position.0 as usize, current_mouse_position.1 as usize,  vx as f64, vy as f64);
                
                pvx += vx.abs();
                pvy += vy.abs();
            }

            last_mouse_position = Some(current_mouse_position);
        } else {
            pvx = 0.0;
            pvy = 0.0;


            last_mouse_position = None;
        }

        texture.update(&image);

        draw_texture_ex(
            texture,
            0.0,
            0.0,
            WHITE,
            DrawTextureParams {
                dest_size: Some(vec2(WIDTH as f32, HEIGHT as f32)),
                ..Default::default()
            },
        );

        let fps_text = format!("FPS: {:.2}", fps);
        draw_text(&fps_text, 10.0, 20.0, 20.0, WHITE);



        next_frame().await;
    }
}