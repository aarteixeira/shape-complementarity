use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use sc_rs::sc::{types::Atom, vector3::Vec3, ScCalculator};

#[pyclass]
pub struct ScResult {
    #[pyo3(get)]
    pub sc: f64,
    #[pyo3(get)]
    pub median_distance: f64,
    #[pyo3(get)]
    pub trimmed_area: f64,
    #[pyo3(get)]
    pub atoms_a: usize,
    #[pyo3(get)]
    pub atoms_b: usize,
}

#[pymethods]
impl ScResult {
    fn __repr__(&self) -> String {
        format!(
            "ScResult(sc={:.4}, median_distance={:.4}, trimmed_area={:.2}, atoms_a={}, atoms_b={})",
            self.sc, self.median_distance, self.trimmed_area, self.atoms_a, self.atoms_b
        )
    }
}

/// Compute Lawrence-Colman Shape Complementarity between two atom groups.
///
/// Mirrors the exact behavior of the sc-rs CLI binary. Atom radii are assigned
/// automatically from atom name + residue name; atoms without a known radius
/// are silently dropped by the upstream library.
///
/// Args:
///     coords_a, atom_names_a, residue_names_a: atoms for molecule A
///     coords_b, atom_names_b, residue_names_b: atoms for molecule B
///     parallel: enable Rayon parallelism inside sc-rs (default True;
///               set False when calling from a ProcessPoolExecutor to avoid
///               oversubscription)
#[pyfunction]
#[pyo3(signature = (coords_a, atom_names_a, residue_names_a, coords_b, atom_names_b, residue_names_b, parallel=true))]
fn compute_sc(
    coords_a: Vec<[f64; 3]>,
    atom_names_a: Vec<String>,
    residue_names_a: Vec<String>,
    coords_b: Vec<[f64; 3]>,
    atom_names_b: Vec<String>,
    residue_names_b: Vec<String>,
    parallel: bool,
) -> PyResult<ScResult> {
    let na = coords_a.len();
    let nb = coords_b.len();

    if na == 0 || nb == 0 {
        return Err(PyValueError::new_err(
            "each atom group must contain at least one atom",
        ));
    }
    if atom_names_a.len() != na || residue_names_a.len() != na {
        return Err(PyValueError::new_err(
            "coords_a, atom_names_a, residue_names_a must all have the same length",
        ));
    }
    if atom_names_b.len() != nb || residue_names_b.len() != nb {
        return Err(PyValueError::new_err(
            "coords_b, atom_names_b, residue_names_b must all have the same length",
        ));
    }

    let mut calc = ScCalculator::new();
    calc.settings_mut().enable_parallel = parallel;

    for i in 0..na {
        let mut a = Atom::new();
        a.coor = Vec3::new(coords_a[i][0], coords_a[i][1], coords_a[i][2]);
        a.atom = atom_names_a[i].clone();
        a.residue = residue_names_a[i].clone();
        calc.add_atom(0, a)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
    }

    for i in 0..nb {
        let mut a = Atom::new();
        a.coor = Vec3::new(coords_b[i][0], coords_b[i][1], coords_b[i][2]);
        a.atom = atom_names_b[i].clone();
        a.residue = residue_names_b[i].clone();
        calc.add_atom(1, a)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
    }

    let results = calc
        .calc()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok(ScResult {
        sc: results.sc,
        median_distance: results.distance,
        trimmed_area: results.area,
        atoms_a: results.surfaces[0].n_atoms,
        atoms_b: results.surfaces[1].n_atoms,
    })
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ScResult>()?;
    m.add_function(wrap_pyfunction!(compute_sc, m)?)?;
    Ok(())
}
