use crate::chromodynamics::Chromodynamics;
use crate::Particle;
use uom::si::f64::{Energy, Force, Length};

enum ForceField {
    Color,
}

impl ForceField {
    fn apply(self, components: Vec<Particle>) -> Vec<Particle> {
        match self {
            ForceField::Color => Chromodynamics::simulate(components),
        }
    }
}
