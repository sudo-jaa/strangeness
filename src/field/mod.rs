use crate::chromodynamics::Chromodynamics;
use crate::system::System;
use crate::Particle;
use uom::si::f64::{Energy, Force, Length, Time};

pub trait Field {
    fn simulate(&self, system: &mut System, step: Time);
}

// pub enum ForceField {
//     Color,
// }
//
// impl ForceField {
//     pub fn apply(self, components: Vec<Particle>) -> System {
//         match self {
//             ForceField::Color => Chromodynamics::simulate(components),
//         }
//     }
// }
