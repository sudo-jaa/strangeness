use crate::field::Field;
use crate::system::System;
use crate::Particle;
use lazy_static::lazy_static;
use ndarray::{Array1, Array2, ArrayBase};
use uom::si::energy::joule;
use uom::si::f64::{Energy, Time};

enum MatrixValue<T> {
    Real(T),
    Imaginary(T),
}

enum SU3Generator {
    λ1,
    λ2,
    λ3,
    λ4,
    λ5,
    λ6,
    λ7,
    λ8,
}

impl SU3Generator {
    fn matrix(&self) -> Array2<i8> {
        match &self {
            SU3Generator::λ1 => {
                todo!()
            }
            SU3Generator::λ2 => {
                todo!()
            }
            SU3Generator::λ3 => {
                todo!()
            }
            SU3Generator::λ4 => {
                todo!()
            }
            SU3Generator::λ5 => {
                todo!()
            }
            SU3Generator::λ6 => {
                todo!()
            }
            SU3Generator::λ7 => {
                todo!()
            }
            SU3Generator::λ8 => {
                todo!()
            }
        }
    }
}

lazy_static! {
    static ref RED_COLUMN_VECTOR: Array1<i8> = Array1::from(vec!(1, 0, 0));
    static ref BLUE_COLUMN_VEC: Array1<i8> = Array1::from(vec!(0, 1, 0));
    static ref GREEN_COLUMN_VEC: Array1<i8> = Array1::from(vec!(0, 0, 1));
    static ref ANTIRED_COLUMN_VECTOR: Array1<i8> = Array1::from(vec!(1, 0, 0));
    static ref ANTIBLUE_COLUMN_VEC: Array1<i8> = Array1::from(vec!(0, 1, 0));
    static ref ANTIGREEN_COLUMN_VEC: Array1<i8> = Array1::from(vec!(0, 0, 1));
}

struct Imaginary(i8);

pub struct Chromodynamics;
impl Field for Chromodynamics {
    fn simulate(&self, system: &mut System, step: Time) {
        // TODO find the matrix representation of the three quark colors and try to intialise one for
        // each colored component in the system. If we have more or less than three we have applied the force to an
        // unstable system.
        println!("Simulating QCD for {:?}", step);

        // TODO
        *system = System {
            particles: vec![],
            free_energy: Default::default(),
        };
    }
}
