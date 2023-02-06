use crate::field::Field;
use crate::Particle;
use uom::si::f64::{Energy, Time};

/// A system containing a number of agent particles and free energy
#[derive(Debug)]
pub struct System {
    pub particles: Vec<Particle>,
    pub free_energy: Energy,
}

impl System {
    /// Allow a system to evolve over a specific time period under the influence of a series
    /// of fields
    pub fn simulate(&mut self, fields: Vec<Box<dyn Field>>, step: Time) {
        fields.iter().for_each(|field| field.simulate(self, step));
    }
}
