use macroquad::prelude::*;
use rayon::prelude::*;

use macroquad::color::{
    hsl_to_rgb
};

const MIN_TEMP: f64 = 0.00000001;
const MAX_TEMP: f64 = 0.005; 

pub struct Fluid {
    dt: f64,
    diff: f64,
    visc: f64,
    width: usize,
    height: usize,
    iter: u32,
}

pub struct FluidGrid {
    width: usize,
    height: usize,
    density: Vec<f64>,
    density_old: Vec<f64>,

    velocity_x: Vec<f64>,
    velocity_y: Vec<f64>,
    
    velocity_x_old: Vec<f64>,
    velocity_y_old: Vec<f64>,    

    
    temperature: Vec<f64>,
    temperature_old: Vec<f64>,

    curl: Vec<f64>,
    color: Vec<Color>
}

impl FluidGrid {
    pub fn new(width: usize, height: usize) -> FluidGrid {
        let size = (width + 2) * (height + 2);
    
        let mut temperature = vec![0.0; size];
        for y in 0..(height + 2) {
            let temp = (height + 1 - y) as f64 / (height + 1) as f64;
            for x in 0..(width + 2) {
                temperature[y * (width + 2) + x] = temp;
            }
        }
    
        FluidGrid {
            width,
            height,
            density: vec![0.0; size],
            velocity_x: vec![0.0; size],
            velocity_y: vec![0.0; size],
    
            temperature,
            temperature_old: vec![0.0; size],
    
            velocity_x_old: vec![0.0; size],
            velocity_y_old: vec![0.0; size],
            density_old: vec![0.0; size],
    
            curl: vec![0.0; size],
            color: vec![Color::new(0.0, 0.0, 0.0, 1.0); size]
        }
    }

    pub fn add_density(&mut self, x: usize, y: usize, amount: f64, color: Color) {
        let idx = self.idx(x, y);
        self.density[idx] += amount;
        self.color[idx] = color;
    }

    pub fn add_velocity(&mut self, x: usize, y: usize, amount_x: f64, amount_y: f64) {
        let idx = self.idx(x, y);
        self.velocity_x[idx] += amount_x;
        self.velocity_y[idx] += amount_y;
    }

    pub fn add_temperature(&mut self, x: usize, y: usize, amount: f64) {
        let idx = self.idx(x, y);

        self.temperature[idx] += amount;
    }

    pub fn apply_buoyancy(&mut self) {
        let buoyancy_coeff = -0.02; 
    
        for j in 1..=self.height {
            for i in 1..=self.width {
                let idx = self.idx(i, j);
                let buoyancy_force = buoyancy_coeff * (self.temperature[idx] - 0.5);
                self.velocity_y[idx] += buoyancy_force;
            }
        }
    }

    pub fn min_max_temperature(&self) -> (f64, f64) {
        let min_temperature = self.temperature.iter().cloned().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let max_temperature = self.temperature.iter().cloned().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    
        (min_temperature, max_temperature)
    }

    pub fn render_densities(&self, image: &mut Image, current_color:&Color, viz_temp: bool) {
        let width = self.width;
        let height = self.height;
        let total_pixels = (width - 2) * (height - 2);


        for idx in 0..total_pixels {
            let x = (idx % (width - 2)) + 1;
            let y = (idx / (width - 2)) + 1;

            let density = self.density[self.idx(x, y)];
            
            let mut color = Color::new(density as f32, 0.1, 0.1, 1.0);

            let mut color = *current_color;
            color.a = density as f32;

            if(viz_temp) {
                let temperature = self.density[self.idx(x,y)].max(0.10).min(0.99);
                color = hsl_to_rgb(temperature as f32, 0.5,0.5);
            }

            image.set_pixel((x - 1) as u32, (y - 1) as u32, color);
        }
    }

    pub fn apply_gravity(&mut self, gravity: f64) {
        let width = self.width;
        let height = self.height;

        self.velocity_y
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, value)| {
                let (i, j) = (idx % (width + 2), idx / (width + 2));
                if i >= 1 && i <= width && j >= 1 && j <= height {
                    *value += gravity;
                }
            });
    }

    fn idx(&self, x: usize, y: usize) -> usize {
        let i = x.min(self.width + 2);
        let j = y.min(self.height + 2);

        return (i + j * (self.width + 2));
    }

}

enum BoundryType {
    X_FIELD,
    Y_FIELD,
    DENSITY_FIELD
}

impl Fluid {
    pub fn new(width: usize, height: usize, dt: f64, diff: f64, visc: f64) -> Fluid {
        let size = (width + 2) * (height + 2);
        Fluid {
            dt,
            diff,
            visc,
            width,
            height,
            iter: 20
        }

    }

    fn idx(&self, x: usize, y: usize) -> usize {
        let i = x.min(self.width + 2);
        let j = y.min(self.height + 2);

        return (i + j * (self.width + 2));
    }


    fn curl(&self, i: usize, j: usize, grid: &FluidGrid) -> f64 {
        (grid.velocity_y[self.idx(i - 1, j)] - grid.velocity_y[self.idx(i + 1, j)]
         + grid.velocity_x[self.idx(i, j + 1)] - grid.velocity_x[self.idx(i, j - 1)])
    }

    pub fn vorticity_confinement(&self, grid: &mut FluidGrid) {
        let mut dw_dx;
        let mut dw_dy;
        let mut length;
        let mut v;
        let vorticity = 5.0;

        for i in 2..(self.width - 1) {
            for j in 2..(self.height - 1) {
                dw_dx = self.curl(i + 0, j - 1, grid).abs() - self.curl(i + 0, j + 1, grid).abs();
                dw_dy = self.curl(i + 1, j, grid).abs() - self.curl(i - 1, j, grid).abs();
    
                length = f64::sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001;
    
                dw_dx = vorticity/length*dw_dx;    
                dw_dy = vorticity/length*dw_dy;
    
                v = self.curl(i, j, grid);
    
                grid.velocity_x[self.idx(i, j)] += dw_dx * v * self.dt;
                grid.velocity_y[self.idx(i, j)] += dw_dy * v * self.dt;
            }
        }
    }

    pub fn step(&self,grid: &mut FluidGrid, vort: bool) {

        self.diffuse(&BoundryType::X_FIELD, &mut grid.velocity_x_old,&grid.velocity_x, self.visc);
        self.diffuse(&BoundryType::Y_FIELD, &mut grid.velocity_y_old, &grid.velocity_y, self.visc);

        self.project(&mut grid.velocity_x_old, &mut grid.velocity_y_old, &mut grid.velocity_x, &mut grid.velocity_y);

        self.advect(&BoundryType::X_FIELD, &mut grid.velocity_x, &grid.velocity_x_old, &grid.velocity_x_old, &grid.velocity_y_old);
        self.advect(&BoundryType::Y_FIELD, &mut grid.velocity_y, &grid.velocity_y_old, &grid.velocity_x_old, &grid.velocity_y_old);

        self.project(&mut grid.velocity_x, &mut grid.velocity_y, &mut grid.velocity_x_old, &mut grid.velocity_y_old);

        self.diffuse(&BoundryType::DENSITY_FIELD, &mut grid.density_old, &grid.density, self.diff);
        self.advect(&BoundryType::DENSITY_FIELD, &mut grid.density, &grid.density_old, &grid.velocity_x, &grid.velocity_y);


        self.fade(&mut grid.density);
        if(vort) {self.vorticity_confinement(grid)};
    }

    fn fade(&self, arr: &mut Vec<f64>) {
        let fade_factor = 1.0/50.0;
        for i in 0..arr.len() {
            arr[i] -= (arr[i] * fade_factor);
        }
    }

    fn project(&self, velocity_x: &mut Vec<f64>, velocity_y: &mut Vec<f64>, p: &mut Vec<f64>, div: &mut Vec<f64>) {
        for j in 1..=self.height {
            for i in 1..=self.width {
                div[self.idx(i, j)] = -0.5 * (velocity_x[self.idx(i + 1, j)] - velocity_x[self.idx(i - 1, j)] + velocity_y[self.idx(i, j + 1)] - velocity_y[self.idx(i, j - 1)]) / self.width as f64;
                p[self.idx(i, j)] = 0.0;
            }
        }
        
        self.set_boundry(&BoundryType::DENSITY_FIELD, div);
        self.set_boundry(&BoundryType::DENSITY_FIELD, p);
        self.lin_solve(&BoundryType::DENSITY_FIELD, p, div, 1.0, 4.0);

        for j in 1..=self.height {
            for i in 1..=self.width {
                velocity_x[self.idx(i, j)] -= 0.5 * self.width as f64 * (p[self.idx(i + 1, j)] - p[self.idx(i - 1, j)]);
                velocity_y[self.idx(i, j)] -= 0.5 * self.width as f64 * (p[self.idx(i, j + 1)] - p[self.idx(i, j - 1)]);
            }
        }

        self.set_boundry(&BoundryType::X_FIELD, velocity_x);
        self.set_boundry(&BoundryType::Y_FIELD, velocity_y);

    }

    fn advect(&self, b: &BoundryType, d: &mut Vec<f64>, d0: &Vec<f64>, velocity_x: &Vec<f64>, velocity_y: &Vec<f64>) {
        let dt0 = self.dt * self.width as f64;
    
        d.par_iter_mut()
            .enumerate()
            .for_each(|(idx, d_value)| {
                let (i, j) = (idx % (self.width + 2), idx / (self.width + 2));
    
                if i >= 1 && i <= self.width && j >= 1 && j <= self.height {
                    let x = (i as f64 - dt0 * velocity_x[idx]).max(0.5).min(self.width as f64 + 0.5);
                    let y = (j as f64 - dt0 * velocity_y[idx]).max(0.5).min(self.height as f64 + 0.5);
    
                    let i0 = x.floor() as usize;
                    let i1 = i0 + 1;
                    let j0 = y.floor() as usize;
                    let j1 = j0 + 1;
    
                    let s1 = x - i0 as f64;
                    let s0 = 1.0 - s1;
                    let t1 = y - j0 as f64;
                    let t0 = 1.0 - t1;
    
                    *d_value = s0 * (t0 * d0[self.idx(i0, j0)] + t1 * d0[self.idx(i0, j1)]) + s1 * (t0 * d0[self.idx(i1, j0)] + t1 * d0[self.idx(i1, j1)]);
                }
            });
    
        self.set_boundry(b, d);
    }

    fn diffuse(&self, b: &BoundryType, x: &mut Vec<f64>, x0: &Vec<f64>, diff: f64) {
        let a = self.dt * diff * (self.width * self.height) as f64;
        self.lin_solve(b, x, x0, a, 1.0 + 4.0*a);    
    }

    // solve linear system using jacobi method
    fn lin_solve(&self, b: &BoundryType, x: &mut Vec<f64>, x0: &Vec<f64>, a: f64, c: f64) {
        let width = self.width;
        let height = self.height;
        let tolerance = 0.00001;
        let n = width * height;
    
        let residual: f64 = x
            .par_iter_mut()
            .enumerate()
            .map(|(idx, x_value)| {
                let (i, j) = (idx % (width + 2), idx / (width + 2));
                if i >= 1 && i <= width && j >= 1 && j <= height {
                    let x_new_value = (x0[idx]
                        + a * (x0[self.idx(i + 1, j)]
                            + x0[self.idx(i - 1, j)]
                            + x0[self.idx(i, j + 1)]
                            + x0[self.idx(i, j - 1)]))
                        / c;
                    let diff = (x_new_value - *x_value).abs();
                    *x_value = x_new_value;
                    diff
                } else {
                    0.0
                }
            })
            .sum();
    
        self.set_boundry(b, x);
    
        if residual / n as f64 > tolerance {
            self.lin_solve(b, x, &x.clone(), a, c);
        }
    }

    fn set_boundry(&self, b: &BoundryType, x: &mut Vec<f64>) {
        for j in 1..=self.height {
            x[self.idx(0,j)] = match b {
                BoundryType::X_FIELD => -x[self.idx(1, j)],
                _ => x[self.idx(1, j)]
            };

            x[self.idx(self.width + 1, j)] = match b {
                BoundryType::X_FIELD => -x[self.idx(self.width, j)],
                _ => x[self.idx(self.width, j)]
            };
        }

        for i in 1..=self.width {
            x[self.idx(i,0)] = match b {
                BoundryType::Y_FIELD => -x[self.idx(i, 1)],
                _ => x[self.idx(i, 1)]
            };

            x[self.idx(i, self.height + 1)] = match b {
                BoundryType::Y_FIELD => -x[self.idx(i, self.height)],
                _ => x[self.idx(i, self.height)]
            };
        }

        x[self.idx(0, 0)] = 0.5 * (x[self.idx(1, 0)] + x[self.idx(0, 1)]);
        x[self.idx(0, self.height + 1)] = 0.5 * (x[self.idx(1, self.height + 1)] + x[self.idx(0, self.height)]);
        x[self.idx(self.width + 1, 0)] = 0.5 * (x[self.idx(self.width, 0)] + x[self.idx(self.width + 1, 1)]);
        x[self.idx(self.width + 1, self.height + 1)] = 0.5 * (x[self.idx(self.width, self.height + 1)] + x[self.idx(self.width + 1, self.height)]);
    }

}