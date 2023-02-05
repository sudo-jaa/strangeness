mod chromodynamics;
mod constants;
mod electrodynamics;
mod field;

use crate::constants::e;
use constants::{h, C, G};
use std::ops::{Add, Neg, Sub};
use uom::num_traits::Pow;
use uom::si::electric_charge::coulomb;
use uom::si::energy::{gigaelectronvolt, joule, megaelectronvolt};
use uom::si::f64::{ElectricCharge, Energy, Force, Length, Mass, Velocity};
use uom::si::force::newton;
use uom::si::mass::kilogram;
use uom::si::velocity::kilometer_per_second;

#[derive(Debug)]
/// Representation of a color charge. This is TBC
enum Color {
    Positive,
    Negative,
}

/// Represents possible combinations of color charges in quarks and gluons. This is just a model,
/// but the internal representation needs to be replaced with SU(3) Glen-mann matrix.
pub type ColorCharge = Option<Vec<Color>>;

#[derive(Debug)]
/// Represents a fractional value
pub struct Fraction(f64);

#[derive(Debug)]
/// Opaque type that represents a spin value as a fraction
pub struct Spin(Fraction);

#[derive(Debug)]
/// Opaque type that represents the charge of a particle as a fraction
struct ParticleCharge(Fraction);

#[derive(Debug)]
/// Represents the mass of a particle.
/// Often in particle physics, the energy and mass of a given particle is used in equivalence.
/// Because of mass-energy equivalence, you can reason about the mass and energy of a particle as a single value
/// that can be computed on the fly from the value you actually have.
pub enum ParticleMass {
    Energy(Energy),
    Mass(Mass),
}

impl ParticleCharge {
    /// Get the electric charge of this particle
    fn to_charge(&self) -> ElectricCharge {
        ElectricCharge::new::<coulomb>(&self.0 .0 * e)
    }

    /// Get the fractional representation of a particle's charge from the real value of an electric charge
    fn from_charge(charge: ElectricCharge) -> ParticleCharge {
        ParticleCharge(Fraction(charge.value / e))
    }
}

#[derive(Debug, Copy, Clone)]
/// Represents possible quantum flavors for a given particle
pub struct QuantumFlavors {
    pub strangeness: i8,
    pub charm: i8,
    pub topness: i8,
    pub bottomness: i8,
}

impl Add for QuantumFlavors {
    type Output = QuantumFlavors;

    fn add(self, rhs: Self) -> Self::Output {
        QuantumFlavors {
            strangeness: self.strangeness + rhs.strangeness,
            charm: self.charm + rhs.charm,
            topness: self.topness + rhs.topness,
            bottomness: self.bottomness + rhs.bottomness,
        }
    }
}

impl Default for QuantumFlavors {
    fn default() -> Self {
        QuantumFlavors {
            strangeness: 0,
            charm: 0,
            topness: 0,
            bottomness: 0,
        }
    }
}

// Vec::with_capacity

impl QuantumFlavors {
    /// Combine any number of quantum flavor systems into a composite quantum flavor system
    fn combine(flavors: Vec<QuantumFlavors>) -> QuantumFlavors {
        flavors
            .iter()
            .fold(QuantumFlavors::default(), |acc, flavor| acc + *flavor)
    }

    /// Flips the state of a quantum flavor system.
    /// This is useful when initialising antiparticles from regular particles
    fn invert(&self) -> QuantumFlavors {
        QuantumFlavors {
            topness: self.topness.neg(),
            charm: self.charm.neg(),
            strangeness: self.strangeness.neg(),
            bottomness: self.bottomness.neg(),
        }
    }
}

/// TODO: output of particle addition
#[derive(Debug)]
struct ParticleCombinationOutput {
    particles: Option<Vec<Particle>>,
    energy: Energy,
}

#[derive(Debug)]
/// A series of properties that describe a fundamental particle
pub struct Particle {
    pub color: ColorCharge,
    pub spin: Spin,
    pub charge: ElectricCharge,
    pub mass: ParticleMass,
    pub flavors: QuantumFlavors,
}

impl Particle {
    fn combine(particles: Vec<Particle>) -> ParticleCombinationOutput {
        todo!()
    }

    /// Produces an anti-particle of a given particle
    fn antiparticle(self) -> Particle {
        Particle::new(
            // needs inverting
            self.color,
            // needs inverting
            self.spin,
            // needs inverting
            ParticleCharge::from_charge(self.charge),
            self.mass,
            self.flavors.invert(),
        )
    }

    /// Initialises a new particle from defaults that represent an up quark
    fn up_quark() -> Particle {
        let mut color: Vec<Color> = Vec::with_capacity(1);
        color.push(Color::Positive);
        Particle::new(
            Some(color),
            Spin(Fraction(1. / 2.)),
            ParticleCharge(Fraction(2. / 3.)),
            ParticleMass::Energy(Energy::new::<megaelectronvolt>(2.2)),
            QuantumFlavors::default(),
        )
    }

    /// Initialises a new particle from defaults that represent a down quark
    fn down_quark() -> Particle {
        let mut color: Vec<Color> = Vec::with_capacity(1);
        color.push(Color::Positive);
        Particle::new(
            Some(color),
            Spin(Fraction(1. / 2.)),
            ParticleCharge(Fraction(-(1. / 3.))),
            ParticleMass::Energy(Energy::new::<megaelectronvolt>(4.7)),
            QuantumFlavors::default(),
        )
    }

    /// Initialises a new particle from defaults that represent a charm quark
    fn charm_quark() -> Particle {
        let mut color: Vec<Color> = Vec::with_capacity(1);
        color.push(Color::Positive);
        Particle::new(
            Some(color),
            Spin(Fraction(1. / 2.)),
            ParticleCharge(Fraction((2. / 3.))),
            ParticleMass::Energy(Energy::new::<gigaelectronvolt>(1.275)),
            QuantumFlavors::default(),
        )
    }

    /// Initialises a new particle from defaults that represent a strange quark
    fn strange_quark() -> Particle {
        let mut color: Vec<Color> = Vec::with_capacity(1);
        color.push(Color::Positive);
        Particle::new(
            Some(color),
            Spin(Fraction(1. / 2.)),
            ParticleCharge(Fraction(-(1. / 3.))),
            ParticleMass::Energy(Energy::new::<megaelectronvolt>(95.)),
            QuantumFlavors::default(),
        )
    }

    /// Initialises a new particle from defaults that represent a top quark
    fn top_quark() -> Particle {
        let mut color: Vec<Color> = Vec::with_capacity(1);
        color.push(Color::Positive);
        Particle::new(
            Some(color),
            Spin(Fraction(1. / 2.)),
            ParticleCharge(Fraction((2. / 3.))),
            ParticleMass::Energy(Energy::new::<gigaelectronvolt>(172.76)),
            QuantumFlavors::default(),
        )
    }

    /// Initialises a new particle from defaults that represent a bottom quark
    fn bottom_quark() -> Particle {
        let mut color: Vec<Color> = Vec::with_capacity(1);
        color.push(Color::Positive);
        Particle::new(
            Some(color),
            Spin(Fraction(1. / 2.)),
            ParticleCharge(Fraction(-(1. / 3.))),
            ParticleMass::Energy(Energy::new::<gigaelectronvolt>(4.18)),
            QuantumFlavors::default(),
        )
    }

    fn electron() -> Particle {
        todo!()
    }
    fn muon() -> Particle {
        todo!()
    }
    fn tau() -> Particle {
        todo!()
    }
    fn electron_neutrino() -> Particle {
        todo!()
    }
    fn muon_neutrino() -> Particle {
        todo!()
    }
    fn tau_neutrino() -> Particle {
        todo!()
    }

    /// Initialise a new particle
    fn new(
        color: ColorCharge,
        spin: Spin,
        charge: ParticleCharge,
        mass: ParticleMass,
        flavors: QuantumFlavors,
    ) -> Particle {
        Particle {
            color,
            spin,
            mass,
            flavors,
            charge: charge.to_charge(),
        }
    }
}

impl Massive for Particle {
    fn get_mass(&self) -> Mass {
        match &self.mass {
            ParticleMass::Energy(energy) => Energetics::calculate_mass(*energy),
            ParticleMass::Mass(mass) => *mass,
        }
    }
}

impl Energetic for Particle {
    fn get_energy(&self) -> Energy {
        match &self.mass {
            ParticleMass::Energy(energy) => *energy,
            ParticleMass::Mass(mass) => Energetics::calculate_energy(*mass),
        }
    }
}

/// Trait indicating entities that have mass
pub trait Massive {
    /// Retrieve the mass of this entity
    fn get_mass(&self) -> Mass;
}

/// Trait indicating entities that have energy
pub trait Energetic {
    /// Retrieve the energy of this entity
    fn get_energy(&self) -> Energy;
}

/// Utilities for working with mass/energy equivalence
pub struct Energetics;
impl Energetics {
    /// Calculate the energy equivalence of mass using E=MC2
    fn calculate_energy(mass: Mass) -> Energy {
        Energy::new::<joule>(mass.value * f64::pow(C, 2))
    }
    /// Calculate the mass equivalence of energy using E=MC2
    fn calculate_mass(energy: Energy) -> Mass {
        Mass::new::<kilogram>(energy.value / f64::pow(C, 2))
    }
}

pub struct Gravitation;

impl Gravitation {
    /// Calculate the gravitational force between two massive objects
    fn calculate_force<T: Massive>(first: T, second: T, distance: Length) -> Force {
        let product_mass = first.get_mass().value * second.get_mass().value;
        let mass_modified_by_g = G * product_mass;
        Force::new::<newton>(mass_modified_by_g / distance.value)
    }
}

#[cfg(test)]
mod tests {
    use crate::Color::Positive;
    use crate::{e, Color, ColorCharge, Fraction, Particle, ParticleCharge, ParticleMass, Spin};
    use approx::relative_eq;
    use uom::si::electric_charge::coulomb;
    use uom::si::electric_charge::Units::franklin;
    use uom::si::energy::megaelectronvolt;
    use uom::si::f64::{ElectricCharge, Energy};

    #[test]
    fn add_fundamental_particles() {
        let up = Particle::up_quark();
        let d1 = Particle::down_quark();
        let d2 = Particle::down_quark();
        let compound = Particle::combine(vec![up, d1, d2]);
        println!("CMP {:?}", compound);
    }

    #[test]
    /// Test to ensure that conversions between fractional and real representations of charges work
    fn fraction_from_real_charge() {
        let charge = ParticleCharge(Fraction(2. / 3.)).to_charge();
        let frac = ParticleCharge::from_charge(charge);
        let charge_2 = frac.to_charge();

        assert_eq!(charge.value, charge_2.value);
    }

    #[test]
    /// Test to ensure that particle will produce a neutral charge in compound
    fn particle_charge_cancellation() {
        let charge_1 = ParticleCharge(Fraction(-(1. / 3.))).to_charge();
        let charge_2 = ParticleCharge(Fraction(-(1. / 3.))).to_charge();
        let charge_3 = ParticleCharge(Fraction(2. / 3.)).to_charge();
        assert_eq!(
            charge_1 + charge_2 + charge_3,
            ElectricCharge::new::<coulomb>(0.)
        );

        let charge_1 = ParticleCharge(Fraction((2. / 3.))).to_charge();
        let charge_2 = ParticleCharge(Fraction((2. / 3.))).to_charge();
        let charge_3 = ParticleCharge(Fraction(-(1. / 3.))).to_charge();
        relative_eq!(
            (charge_1 + charge_2 + charge_3).value,
            ElectricCharge::new::<coulomb>(1. * e).value
        );
    }
}
