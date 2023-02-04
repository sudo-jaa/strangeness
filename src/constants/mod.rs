/// Gravitational constant
pub const G: f64 = 6.67408e-11;
// Speed of Light
pub const C: f64 = 299792458.0;
// Planck Constant
pub const h: f64 = 6.62607015e-34;
// Elementary charge
pub const e: f64 = 1.602176634e-19;

#[cfg(test)]
mod tests {
    use approx::relative_eq;
    use uom::si::f64::{Length, Mass, Energy, Force};
    use uom::si::length::meter;
    use uom::si::mass::kilogram;
    use uom::si::energy::joule;
    use uom::si::force::newton;
    use crate::{Gravitation, Massive, Energetics};
    use super::*;


    #[test]
    fn energy_mass_equivalence_is_correct() {
        struct TestObject {
            mass: Mass
        }

        impl Massive for TestObject {
            fn get_mass(&self) -> Mass {
                self.mass
            }
        }

        let object = TestObject{mass: Mass::new::<kilogram>(0.000000000001)};
        relative_eq!(Energetics::calculate_energy(object.mass).value, Energy::new::<joule>(89875.5).value);

        let object = TestObject{mass: Mass::new::<kilogram>(1.)};
        relative_eq!(Energetics::calculate_energy(object.mass).value, Energy::new::<joule>(89875517873681764.).value);
    }

    #[test]
    fn gravitational_constant_calculations_are_correct() {
        struct TestObject {
            mass: Mass
        }

        impl Massive for TestObject {
            fn get_mass(&self) -> Mass {
                self.mass
            }
        }

        let gravitational_force = (G * 1. * 1.) / 1.0;
        assert_eq!(gravitational_force, G);

        let first = TestObject { mass: Mass::new::<kilogram>(1.0) };
        let second = TestObject { mass: Mass::new::<kilogram>(1.0) };
        assert_eq!(Gravitation::calculate_force(first, second, Length::new::<meter>(1.0)).value, gravitational_force);
    }

    #[test]
    fn gravitational_calculations_are_correct() {
        struct TestObject {
            mass: Mass
        }

        impl Massive for TestObject {
            fn get_mass(&self) -> Mass {
                self.mass
            }
        }

        let first = TestObject { mass: Mass::new::<kilogram>(10.0) };
        let second = TestObject { mass: Mass::new::<kilogram>(15.0) };
        relative_eq!(Gravitation::calculate_force(first, second, Length::new::<meter>(1.0)).value, Force::new::<newton>(0.00000001001145).value);

        let first = TestObject { mass: Mass::new::<kilogram>(1210.1) };
        let second = TestObject { mass: Mass::new::<kilogram>(125.4) };
        relative_eq!(Gravitation::calculate_force(first, second, Length::new::<meter>(12.6)).value, Force::new::<newton>(0.0000000637945).value);

        let first = TestObject { mass: Mass::new::<kilogram>(1989000000000000000000000000000.) };
        let second = TestObject { mass: Mass::new::<kilogram>(125.4) };
        relative_eq!(Gravitation::calculate_force(first, second, Length::new::<meter>(695700000.)).value, Force::new::<newton>(34394.9).value);

    }
}

