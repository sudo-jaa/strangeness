use crate::Particle;
use lazy_static::lazy_static;
use ndarray::{Array2, ArrayBase};

enum MatrixValue<T> {
    Real(T),
    Imaginary(T),
}

lazy_static! {
    static ref VALID_GELLMANN_MATRICES: [Array2<MatrixValue>; 8] = [
        Array2::from(vec!(
            [
                MatrixValue::Real(0),
                MatrixValue::Real(1),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(1),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ]
        )),
        Array2::from(vec!(
            [
                MatrixValue::Real(0),
                MatrixValue::Imaginary(-1),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Imaginary(1),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ]
        )),
        Array2::from(vec!(
            [
                MatrixValue::Real(1),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(-1),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ]
        )),
        Array2::from(vec!(
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(1)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(1),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ]
        )),
        Array2::from(vec!(
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Imaginary(-1)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Imaginary(1),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ]
        )),
        Array2::from(vec!(
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(1)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(1),
                MatrixValue::Real(0)
            ]
        )),
        Array2::from(vec!(
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Imaginary(-1)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Imaginary(1),
                MatrixValue::Real(0)
            ]
        )),
        Array2::from(vec!(
            [
                MatrixValue::Real(1),
                MatrixValue::Real(0),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(1),
                MatrixValue::Real(0)
            ],
            [
                MatrixValue::Real(0),
                MatrixValue::Real(0),
                MatrixValue::Real(-2)
            ]
        )),
    ];
}

struct Imaginary(i8);

pub struct Chromodynamics;
impl Chromodynamics {
    pub fn simulate(components: Vec<Particle>) -> Vec<Particle> {
        // TODO find the matrix representation of the three quark colors and try to intialise one for
        // each component in the system. If we have more or less than three we have applied the force to an
        // unstable system.
        let gellmann_matrices = components
            .iter()
            .map(|particle| {
                if particle.color.is_none() {
                    None
                } else {
                    todo!()
                }
            })
            .collect::<Vec<Option<Array2<MatrixValue<i8>>>>>();
        todo!()
    }
}
