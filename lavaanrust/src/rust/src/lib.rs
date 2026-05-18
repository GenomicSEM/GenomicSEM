use extendr_api::prelude::*;
use nalgebra::{DMatrix, DVector, SymmetricEigen};

fn from_col_major(values: &[f64], nrows: usize, ncols: usize) -> DMatrix<f64> {
    DMatrix::from_column_slice(nrows, ncols, values)
}

fn to_col_major(matrix: &DMatrix<f64>) -> Vec<f64> {
    matrix.as_slice().to_vec()
}

fn vech(matrix: &DMatrix<f64>) -> DVector<f64> {
    let n = matrix.nrows();
    let mut values = Vec::with_capacity(n * (n + 1) / 2);

    for col in 0..n {
        for row in col..n {
            values.push(matrix[(row, col)]);
        }
    }

    DVector::from_vec(values)
}

enum DwlsWeights {
    Diagonal(DVector<f64>),
    Dense(DMatrix<f64>),
}

impl DwlsWeights {
    fn from_col_major(values: &[f64], n_stats: usize) -> Self {
        let matrix = from_col_major(values, n_stats, n_stats);
        if (0..n_stats).all(|col| (0..n_stats).all(|row| row == col || matrix[(row, col)] == 0.0)) {
            Self::Diagonal(matrix.diagonal())
        } else {
            Self::Dense(matrix)
        }
    }

    fn apply_vector(&self, vector: &DVector<f64>) -> DVector<f64> {
        match self {
            Self::Diagonal(diagonal) => diagonal.component_mul(vector),
            Self::Dense(matrix) => matrix * vector,
        }
    }

    fn apply_matrix(&self, matrix: &DMatrix<f64>) -> DMatrix<f64> {
        match self {
            Self::Diagonal(diagonal) => {
                let mut weighted = matrix.clone();
                for row in 0..weighted.nrows() {
                    weighted.row_mut(row).scale_mut(diagonal[row]);
                }
                weighted
            }
            Self::Dense(weights) => weights * matrix,
        }
    }

    fn objective(&self, residual: &DVector<f64>) -> f64 {
        0.5 * residual.dot(&self.apply_vector(residual))
    }

    fn normal_equations(
        &self,
        delta: &DMatrix<f64>,
        residual: &DVector<f64>,
    ) -> (DVector<f64>, DMatrix<f64>) {
        match self {
            Self::Diagonal(diagonal) => diagonal_normal_equations(delta, diagonal, Some(residual))
                .unwrap_or_else(|| dense_normal_equations(self, delta, residual)),
            Self::Dense(_) => dense_normal_equations(self, delta, residual),
        }
    }

    fn weighted_crossprod(&self, delta: &DMatrix<f64>) -> DMatrix<f64> {
        match self {
            Self::Diagonal(diagonal) => diagonal_normal_equations(delta, diagonal, None)
                .map(|(_, hessian)| hessian)
                .unwrap_or_else(|| {
                    let weighted_delta = self.apply_matrix(delta);
                    delta.transpose() * &weighted_delta
                }),
            Self::Dense(_) => {
                let weighted_delta = self.apply_matrix(delta);
                delta.transpose() * &weighted_delta
            }
        }
    }
}

fn dense_normal_equations(
    weights: &DwlsWeights,
    delta: &DMatrix<f64>,
    residual: &DVector<f64>,
) -> (DVector<f64>, DMatrix<f64>) {
    let weighted_residual = weights.apply_vector(residual);
    let weighted_delta = weights.apply_matrix(delta);
    (
        delta.transpose() * &weighted_residual,
        delta.transpose() * &weighted_delta,
    )
}

fn diagonal_normal_equations(
    delta: &DMatrix<f64>,
    diagonal: &DVector<f64>,
    residual: Option<&DVector<f64>>,
) -> Option<(DVector<f64>, DMatrix<f64>)> {
    let n_rows = delta.nrows();
    let n_cols = delta.ncols();
    let mut row_entries = Vec::with_capacity(n_rows);
    let mut sparse_pair_work = 0usize;

    for row in 0..n_rows {
        let mut entries = Vec::new();
        for col in 0..n_cols {
            let value = delta[(row, col)];
            if value != 0.0 {
                entries.push((col, value));
            }
        }

        sparse_pair_work = sparse_pair_work.saturating_add(entries.len() * (entries.len() + 1) / 2);
        row_entries.push(entries);
    }

    let dense_pair_work = n_rows.saturating_mul(n_cols).saturating_mul(n_cols);
    if sparse_pair_work.saturating_mul(8) >= dense_pair_work {
        return None;
    }

    let mut gradient = DVector::<f64>::zeros(n_cols);
    let mut hessian = DMatrix::<f64>::zeros(n_cols, n_cols);

    for row in 0..n_rows {
        let row_weight = diagonal[row];
        let entries = &row_entries[row];

        if let Some(residual) = residual {
            let weighted_residual = row_weight * residual[row];
            for &(col, value) in entries {
                gradient[col] += value * weighted_residual;
            }
        }

        for (left_pos, &(left_col, left_value)) in entries.iter().enumerate() {
            for &(right_col, right_value) in &entries[..=left_pos] {
                let contribution = row_weight * left_value * right_value;
                hessian[(left_col, right_col)] += contribution;
                if left_col != right_col {
                    hessian[(right_col, left_col)] += contribution;
                }
            }
        }
    }

    Some((gradient, hessian))
}

#[cfg(test)]
mod dwls_weight_tests {
    use super::*;

    #[test]
    fn diagonal_weights_scale_without_dense_multiplication() {
        let weights = DwlsWeights::from_col_major(
            &[
                2.0, 0.0, 0.0, //
                0.0, 3.0, 0.0, //
                0.0, 0.0, 4.0,
            ],
            3,
        );
        let vector = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let matrix = DMatrix::from_row_slice(3, 2, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);

        assert!(matches!(&weights, DwlsWeights::Diagonal(_)));
        assert_eq!(
            weights.apply_vector(&vector),
            DVector::from_vec(vec![2.0, 6.0, 12.0])
        );
        assert_eq!(
            weights.apply_matrix(&matrix),
            DMatrix::from_row_slice(3, 2, &[2.0, 4.0, 9.0, 12.0, 20.0, 24.0])
        );
    }

    #[test]
    fn dense_weights_preserve_general_matrix_products() {
        let values = vec![
            2.0, 0.5, 0.0, //
            0.5, 3.0, 0.25, //
            0.0, 0.25, 4.0,
        ];
        let weights = DwlsWeights::from_col_major(&values, 3);
        let dense = from_col_major(&values, 3, 3);
        let vector = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let matrix = DMatrix::from_row_slice(3, 2, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);

        assert!(matches!(&weights, DwlsWeights::Dense(_)));
        assert_eq!(weights.apply_vector(&vector), &dense * &vector);
        assert_eq!(weights.apply_matrix(&matrix), &dense * &matrix);
    }

    #[test]
    fn sparse_diagonal_normal_equations_match_dense_products() {
        let weights = DwlsWeights::Diagonal(DVector::from_vec(vec![2.0, 3.0, 4.0, 5.0]));
        let residual = DVector::from_vec(vec![1.0, -2.0, 0.5, 3.0]);
        let delta = DMatrix::from_row_slice(
            4,
            4,
            &[
                1.0, 0.0, 0.0, 0.0, //
                0.0, 2.0, 0.0, 0.0, //
                0.0, 0.0, 3.0, 0.0, //
                0.0, 0.0, 0.0, 4.0,
            ],
        );
        let dense = dense_normal_equations(&weights, &delta, &residual);
        let sparse = weights.normal_equations(&delta, &residual);

        assert_eq!(sparse.0, dense.0);
        assert_eq!(sparse.1, dense.1);
        assert_eq!(weights.weighted_crossprod(&delta), dense.1);
    }
}

#[cfg(test)]
mod ram_surface_tests {
    use super::*;

    #[test]
    fn sparse_rank_one_vech_matches_dense_symmetric_outer() {
        let left = DMatrix::from_column_slice(4, 1, &[1.0, 0.0, 2.0, 0.0]);
        let right = DMatrix::from_column_slice(4, 1, &[0.0, 3.0, 4.0, 0.0]);
        let mut sparse = DVector::zeros(10);
        let left_support = column_supports(&left);
        let right_support = column_supports(&right);

        add_symmetric_outer_vech(
            &mut sparse.as_view_mut(),
            &left,
            0,
            &left_support[0],
            &right,
            0,
            &right_support[0],
            4,
        );

        let dense = &left * right.transpose() + &right * left.transpose();
        assert_eq!(sparse, vech(&dense));
    }

    #[test]
    fn sparse_self_outer_vech_matches_dense_outer() {
        let matrix = DMatrix::from_column_slice(4, 1, &[1.0, 0.0, 2.0, 0.0]);
        let mut sparse = DVector::zeros(10);
        let support = column_supports(&matrix);

        add_self_outer_vech(&mut sparse.as_view_mut(), &matrix, 0, &support[0], 4);

        let dense = &matrix * matrix.transpose();
        assert_eq!(sparse, vech(&dense));
    }

    #[test]
    fn damped_normal_equation_solver_matches_dense_solution() {
        let matrix = DMatrix::from_row_slice(3, 3, &[4.0, 1.0, 0.0, 1.0, 3.0, 0.5, 0.0, 0.5, 2.0]);
        let rhs = DVector::from_vec(vec![1.0, 2.0, 3.0]);

        let expected = matrix.clone().lu().solve(&rhs).unwrap();
        let actual = solve_damped_normal_equations(&matrix, &rhs).unwrap();

        assert!((&actual - expected).amax() < 1e-12);
    }
}

fn implied_one_factor(loadings: &DVector<f64>, residuals: &DVector<f64>) -> DMatrix<f64> {
    let mut implied = loadings * loadings.transpose();

    for idx in 0..residuals.len() {
        implied[(idx, idx)] += residuals[idx];
    }

    implied
}

fn implied_commonfactor_gwas(
    loadings: &DVector<f64>,
    gamma: f64,
    residuals: &DVector<f64>,
    psi: f64,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let mut implied = DMatrix::<f64>::zeros(k + 1, k + 1);
    let factor_total = gamma * gamma * phi + psi;

    implied[(0, 0)] = phi;

    for idx in 0..k {
        let cov = phi * loadings[idx] * gamma;
        implied[(idx + 1, 0)] = cov;
        implied[(0, idx + 1)] = cov;
    }

    for col in 0..k {
        for row in col..k {
            let mut cov = loadings[row] * loadings[col] * factor_total;
            if row == col {
                cov += residuals[row];
            }

            implied[(row + 1, col + 1)] = cov;
            implied[(col + 1, row + 1)] = cov;
        }
    }

    implied
}

fn delta_commonfactor_gwas(
    loadings: &DVector<f64>,
    gamma: f64,
    psi: f64,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let n_stats = (k + 1) * (k + 2) / 2;
    let n_params = 2 * k + 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, n_params);
    let factor_total = gamma * gamma * phi + psi;
    let gamma_col = k - 1;
    let theta_start = gamma_col + 1;
    let psi_col = theta_start + k;
    let phi_col = psi_col + 1;
    let mut stat_row = 0;

    for col in 0..=k {
        for row in col..=k {
            if col == 0 && row == 0 {
                delta[(stat_row, phi_col)] = 1.0;
            } else if col == 0 || row == 0 {
                let trait_idx = row.max(col) - 1;

                if trait_idx > 0 {
                    delta[(stat_row, trait_idx - 1)] = phi * gamma;
                }

                delta[(stat_row, gamma_col)] = phi * loadings[trait_idx];
                delta[(stat_row, phi_col)] = loadings[trait_idx] * gamma;
            } else {
                let trait_row = row - 1;
                let trait_col = col - 1;

                if trait_row > 0 {
                    delta[(stat_row, trait_row - 1)] += loadings[trait_col] * factor_total;
                }

                if trait_col > 0 {
                    delta[(stat_row, trait_col - 1)] += loadings[trait_row] * factor_total;
                }

                delta[(stat_row, gamma_col)] =
                    loadings[trait_row] * loadings[trait_col] * 2.0 * gamma * phi;
                delta[(stat_row, psi_col)] = loadings[trait_row] * loadings[trait_col];
                delta[(stat_row, phi_col)] =
                    loadings[trait_row] * loadings[trait_col] * gamma * gamma;

                if trait_row == trait_col {
                    delta[(stat_row, theta_start + trait_row)] = 1.0;
                }
            }

            stat_row += 1;
        }
    }

    delta
}

fn implied_commonfactor_gwas_q(
    loadings: &DVector<f64>,
    gamma: f64,
    direct: &DVector<f64>,
    residuals: &DVector<f64>,
    psi: f64,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let mut implied = DMatrix::<f64>::zeros(k + 1, k + 1);
    let beta = loadings * gamma + direct;

    implied[(0, 0)] = phi;

    for idx in 0..k {
        let cov = phi * beta[idx];
        implied[(idx + 1, 0)] = cov;
        implied[(0, idx + 1)] = cov;
    }

    for col in 0..k {
        for row in col..k {
            let mut cov = loadings[row] * loadings[col] * psi + phi * beta[row] * beta[col];
            if row == col {
                cov += residuals[row];
            }

            implied[(row + 1, col + 1)] = cov;
            implied[(col + 1, row + 1)] = cov;
        }
    }

    implied
}

fn delta_commonfactor_gwas_q(
    loadings: &DVector<f64>,
    gamma: f64,
    direct: &DVector<f64>,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let n_stats = (k + 1) * (k + 2) / 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, 2 * k);
    let beta = loadings * gamma + direct;
    let mut stat_row = 0;

    for col in 0..=k {
        for row in col..=k {
            if col == 0 && row == 0 {
                // The predictor variance is fixed in the Q-model refit.
            } else if col == 0 || row == 0 {
                let trait_idx = row.max(col) - 1;
                delta[(stat_row, trait_idx)] = phi;
            } else {
                let trait_row = row - 1;
                let trait_col = col - 1;

                delta[(stat_row, trait_row)] += phi * beta[trait_col];
                delta[(stat_row, trait_col)] += phi * beta[trait_row];

                if trait_row == trait_col {
                    delta[(stat_row, k + trait_row)] = 1.0;
                }
            }

            stat_row += 1;
        }
    }

    delta
}

fn delta_user_gwas_fixed_measurement(
    loadings: &DVector<f64>,
    gamma: f64,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let n_stats = (k + 1) * (k + 2) / 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, k + 3);
    let psi_col = k;
    let gamma_col = k + 1;
    let phi_col = k + 2;
    let mut stat_row = 0;

    for col in 0..=k {
        for row in col..=k {
            if col == 0 && row == 0 {
                delta[(stat_row, phi_col)] = 1.0;
            } else if col == 0 || row == 0 {
                let trait_idx = row.max(col) - 1;
                delta[(stat_row, gamma_col)] = phi * loadings[trait_idx];
                delta[(stat_row, phi_col)] = loadings[trait_idx] * gamma;
            } else {
                let trait_row = row - 1;
                let trait_col = col - 1;

                if trait_row == trait_col {
                    delta[(stat_row, trait_row)] = 1.0;
                }

                delta[(stat_row, psi_col)] = loadings[trait_row] * loadings[trait_col];
                delta[(stat_row, gamma_col)] =
                    loadings[trait_row] * loadings[trait_col] * 2.0 * gamma * phi;
                delta[(stat_row, phi_col)] =
                    loadings[trait_row] * loadings[trait_col] * gamma * gamma;
            }

            stat_row += 1;
        }
    }

    delta
}

fn delta_one_factor(loadings: &DVector<f64>) -> DMatrix<f64> {
    let n = loadings.len();
    let n_stats = n * (n + 1) / 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, 2 * n);
    let mut row = 0;

    for col in 0..n {
        for obs_row in col..n {
            if obs_row == col {
                delta[(row, col)] = 2.0 * loadings[col];
                delta[(row, n + col)] = 1.0;
            } else {
                delta[(row, col)] = loadings[obs_row];
                delta[(row, obs_row)] = loadings[col];
            }

            row += 1;
        }
    }

    delta
}

fn initial_one_factor(sample_cov: &DMatrix<f64>) -> (DVector<f64>, DVector<f64>) {
    let eigen = SymmetricEigen::new(sample_cov.clone());
    let mut best_idx = 0;

    for idx in 1..eigen.eigenvalues.len() {
        if eigen.eigenvalues[idx] > eigen.eigenvalues[best_idx] {
            best_idx = idx;
        }
    }

    let scale = eigen.eigenvalues[best_idx].max(1e-8).sqrt();
    let mut loadings = eigen.eigenvectors.column(best_idx).into_owned() * scale;

    if loadings.sum() < 0.0 {
        loadings *= -1.0;
    }

    let mut residuals = DVector::<f64>::zeros(sample_cov.nrows());
    for idx in 0..sample_cov.nrows() {
        residuals[idx] = (sample_cov[(idx, idx)] - loadings[idx] * loadings[idx]).max(1e-8);
    }

    (loadings, residuals)
}

fn initial_commonfactor_gwas(
    sample_cov: &DMatrix<f64>,
) -> (DVector<f64>, f64, DVector<f64>, f64, f64) {
    let k = sample_cov.nrows() - 1;
    let trait_cov = sample_cov.view((1, 1), (k, k)).into_owned();
    let (std_loadings, residuals) = initial_one_factor(&trait_cov);
    let marker = if std_loadings[0].abs() < 1e-8 {
        1.0
    } else {
        std_loadings[0]
    };
    let mut loadings = std_loadings / marker;
    loadings[0] = 1.0;
    let phi = sample_cov[(0, 0)].max(1e-8);
    let gamma = sample_cov[(1, 0)] / phi;
    let factor_total = marker * marker;
    let psi = (factor_total - gamma * gamma * phi).max(1e-8);

    (loadings, gamma, residuals, psi, phi)
}

fn srmr(sample_cov: &DMatrix<f64>, implied: &DMatrix<f64>) -> f64 {
    let observed_sd: Vec<f64> = (0..sample_cov.nrows())
        .map(|idx| sample_cov[(idx, idx)].abs().sqrt())
        .collect();
    let mut sum_sq = 0.0;
    let mut count = 0usize;

    for col in 0..sample_cov.ncols() {
        for row in col..sample_cov.nrows() {
            let denom = observed_sd[row] * observed_sd[col];
            if denom > 0.0 {
                let resid = (sample_cov[(row, col)] - implied[(row, col)]) / denom;
                sum_sq += resid * resid;
                count += 1;
            }
        }
    }

    if count == 0 {
        0.0
    } else {
        (sum_sq / count as f64).sqrt()
    }
}

fn from_vech(values: &DVector<f64>, n: usize) -> DMatrix<f64> {
    let mut matrix = DMatrix::<f64>::zeros(n, n);
    let mut idx = 0;

    for col in 0..n {
        for row in col..n {
            matrix[(row, col)] = values[idx];
            matrix[(col, row)] = values[idx];
            idx += 1;
        }
    }

    matrix
}

fn select_rows(matrix: &DMatrix<f64>, row_indices: &[usize]) -> DMatrix<f64> {
    let mut selected = DMatrix::<f64>::zeros(row_indices.len(), matrix.ncols());

    for (out_row, source_row) in row_indices.iter().enumerate() {
        for col in 0..matrix.ncols() {
            selected[(out_row, col)] = matrix[(*source_row, col)];
        }
    }

    selected
}

fn implied_ram(
    directed: &DMatrix<f64>,
    covariance: &DMatrix<f64>,
    observed_indices: &[usize],
) -> std::result::Result<(DMatrix<f64>, DMatrix<f64>, DMatrix<f64>), Error> {
    let identity = DMatrix::<f64>::identity(directed.nrows(), directed.ncols());
    let Some(inverse) = (identity - directed).try_inverse() else {
        return Err(Error::Other(
            "RAM solve failed because I - A is singular".into(),
        ));
    };
    let observed_inverse = select_rows(&inverse, observed_indices);
    let full_to_observed = &inverse * covariance * observed_inverse.transpose();
    let implied = select_rows(&full_to_observed, observed_indices);
    let implied_by_variable = full_to_observed.transpose();

    Ok((implied, implied_by_variable, observed_inverse))
}

fn implied_ram_observed(
    directed: &DMatrix<f64>,
    covariance: &DMatrix<f64>,
    observed_indices: &[usize],
) -> std::result::Result<DMatrix<f64>, Error> {
    let identity = DMatrix::<f64>::identity(directed.nrows(), directed.ncols());
    let Some(inverse) = (identity - directed).try_inverse() else {
        return Err(Error::Other(
            "RAM solve failed because I - A is singular".into(),
        ));
    };
    let mut observed_inverse = DMatrix::<f64>::zeros(observed_indices.len(), inverse.ncols());

    for (out_row, source_row) in observed_indices.iter().enumerate() {
        for col in 0..inverse.ncols() {
            observed_inverse[(out_row, col)] = inverse[(*source_row, col)];
        }
    }

    Ok(&observed_inverse * covariance * observed_inverse.transpose())
}

fn validate_ram_indices(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    fixed_values: &[f64],
    observed_index: &[i32],
    n_variables: usize,
    n_free: usize,
) -> std::result::Result<Vec<usize>, Error> {
    let n_rows = lhs_index.len();
    if rhs_index.len() != n_rows
        || op_code.len() != n_rows
        || free_index.len() != n_rows
        || fixed_values.len() != n_rows
    {
        return Err(Error::Other(
            "RAM row vectors must all have the same length".into(),
        ));
    }

    if observed_index.is_empty() {
        return Err(Error::Other("observed_index must be non-empty".into()));
    }

    for row_idx in 0..n_rows {
        let lhs = lhs_index[row_idx];
        let rhs = rhs_index[row_idx];
        if lhs <= 0 || rhs <= 0 || lhs as usize > n_variables || rhs as usize > n_variables {
            return Err(Error::Other("RAM row index is out of bounds".into()));
        }

        let free = free_index[row_idx];
        if free < 0 || free as usize > n_free {
            return Err(Error::Other("free_index is out of bounds".into()));
        }

        if !matches!(op_code[row_idx], 1..=3) {
            return Err(Error::Other("unsupported RAM operator code".into()));
        }
    }

    observed_index
        .iter()
        .map(|index| {
            if *index <= 0 || *index as usize > n_variables {
                Err(Error::Other("observed_index is out of bounds".into()))
            } else {
                Ok((*index as usize) - 1)
            }
        })
        .collect::<std::result::Result<Vec<_>, _>>()
}

fn validate_free_row_groups(
    free_row_offsets: &[i32],
    free_row_indices: &[i32],
    free_index: &[i32],
    n_free: usize,
) -> std::result::Result<(), Error> {
    if free_row_offsets.len() != n_free + 1 {
        return Err(Error::Other(
            "free_row_offsets must have length n_free + 1".into(),
        ));
    }

    if free_row_offsets.first().copied() != Some(0)
        || free_row_offsets.last().copied() != Some(free_row_indices.len() as i32)
    {
        return Err(Error::Other(
            "free_row_offsets must span free_row_indices".into(),
        ));
    }

    for window in free_row_offsets.windows(2) {
        if window[0] > window[1] {
            return Err(Error::Other(
                "free_row_offsets must be non-decreasing".into(),
            ));
        }
    }

    for (free_position, window) in free_row_offsets.windows(2).enumerate() {
        for flat_idx in window[0] as usize..window[1] as usize {
            let row_idx = free_row_indices[flat_idx];
            if row_idx <= 0 || row_idx as usize > free_index.len() {
                return Err(Error::Other("free_row_indices is out of bounds".into()));
            }

            if free_index[row_idx as usize - 1] != free_position as i32 + 1 {
                return Err(Error::Other(
                    "free_row_indices does not match free_index".into(),
                ));
            }
        }
    }

    Ok(())
}

struct RamDwlsPlan {
    lhs_index: Vec<i32>,
    rhs_index: Vec<i32>,
    op_code: Vec<i32>,
    free_index: Vec<i32>,
    fixed_values: Vec<f64>,
    observed_indices: Vec<usize>,
    free_row_offsets: Vec<i32>,
    free_row_indices: Vec<i32>,
    lower_bounds: Vec<f64>,
    upper_bounds: Vec<f64>,
    n_variables: usize,
    n_observed: usize,
    n_stats: usize,
}

impl RamDwlsPlan {
    fn new(
        lhs_index: Vec<i32>,
        rhs_index: Vec<i32>,
        op_code: Vec<i32>,
        free_index: Vec<i32>,
        fixed_values: Vec<f64>,
        observed_index: Vec<i32>,
        free_row_offsets: Vec<i32>,
        free_row_indices: Vec<i32>,
        lower_bounds: Vec<f64>,
        upper_bounds: Vec<f64>,
        n_variables: usize,
    ) -> std::result::Result<Self, Error> {
        if lower_bounds.len() != upper_bounds.len() {
            return Err(Error::Other(
                "lower_bounds and upper_bounds must have the same size".into(),
            ));
        }

        if lower_bounds
            .iter()
            .zip(upper_bounds.iter())
            .any(|(lower, upper)| lower > upper)
        {
            return Err(Error::Other("lower_bounds exceed upper_bounds".into()));
        }

        let observed_indices = validate_ram_indices(
            &lhs_index,
            &rhs_index,
            &op_code,
            &free_index,
            &fixed_values,
            &observed_index,
            n_variables,
            lower_bounds.len(),
        )?;
        validate_free_row_groups(
            &free_row_offsets,
            &free_row_indices,
            &free_index,
            lower_bounds.len(),
        )?;

        let n_observed = observed_indices.len();
        let n_stats = n_observed * (n_observed + 1) / 2;

        Ok(Self {
            lhs_index,
            rhs_index,
            op_code,
            free_index,
            fixed_values,
            observed_indices,
            free_row_offsets,
            free_row_indices,
            lower_bounds,
            upper_bounds,
            n_variables,
            n_observed,
            n_stats,
        })
    }

    fn n_free(&self) -> usize {
        self.lower_bounds.len()
    }
}

fn ram_matrices_from_rows(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    fixed_values: &[f64],
    free_values: &[f64],
    n_variables: usize,
) -> (DMatrix<f64>, DMatrix<f64>) {
    let mut directed = DMatrix::<f64>::zeros(n_variables, n_variables);
    let mut covariance = DMatrix::<f64>::zeros(n_variables, n_variables);

    for row_idx in 0..lhs_index.len() {
        let lhs = lhs_index[row_idx] as usize - 1;
        let rhs = rhs_index[row_idx] as usize - 1;
        let free = free_index[row_idx];
        let value = if free > 0 {
            free_values[free as usize - 1]
        } else {
            fixed_values[row_idx]
        };

        match op_code[row_idx] {
            1 => directed[(rhs, lhs)] = value,
            2 => directed[(lhs, rhs)] = value,
            3 => {
                covariance[(lhs, rhs)] = value;
                covariance[(rhs, lhs)] = value;
            }
            _ => unreachable!(),
        }
    }

    (directed, covariance)
}

fn ram_delta_from_rows(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_row_offsets: &[i32],
    free_row_indices: &[i32],
    implied_by_variable: &DMatrix<f64>,
    observed_inverse: &DMatrix<f64>,
    n_free: usize,
) -> DMatrix<f64> {
    let n_observed = observed_inverse.nrows();
    let n_stats = n_observed * (n_observed + 1) / 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, n_free);
    let implied_supports = column_supports(implied_by_variable);
    let inverse_supports = column_supports(observed_inverse);

    for free_position in 0..n_free {
        let start = free_row_offsets[free_position] as usize;
        let end = free_row_offsets[free_position + 1] as usize;
        let mut column = delta.column_mut(free_position);

        for flat_idx in start..end {
            let row_idx = free_row_indices[flat_idx] as usize - 1;
            let lhs = lhs_index[row_idx] as usize - 1;
            let rhs = rhs_index[row_idx] as usize - 1;

            match op_code[row_idx] {
                1 => {
                    add_symmetric_outer_vech(
                        &mut column,
                        observed_inverse,
                        rhs,
                        &inverse_supports[rhs],
                        implied_by_variable,
                        lhs,
                        &implied_supports[lhs],
                        n_observed,
                    );
                }
                2 => {
                    add_symmetric_outer_vech(
                        &mut column,
                        observed_inverse,
                        lhs,
                        &inverse_supports[lhs],
                        implied_by_variable,
                        rhs,
                        &implied_supports[rhs],
                        n_observed,
                    );
                }
                3 if lhs == rhs => {
                    add_self_outer_vech(
                        &mut column,
                        observed_inverse,
                        lhs,
                        &inverse_supports[lhs],
                        n_observed,
                    );
                }
                3 => {
                    add_symmetric_outer_vech(
                        &mut column,
                        observed_inverse,
                        lhs,
                        &inverse_supports[lhs],
                        observed_inverse,
                        rhs,
                        &inverse_supports[rhs],
                        n_observed,
                    );
                }
                _ => unreachable!(),
            }
        }
    }

    delta
}

fn column_supports(matrix: &DMatrix<f64>) -> Vec<Vec<usize>> {
    (0..matrix.ncols())
        .map(|col| {
            (0..matrix.nrows())
                .filter(|row| matrix[(*row, col)] != 0.0)
                .collect::<Vec<_>>()
        })
        .collect()
}

fn vech_index(row: usize, col: usize, n_observed: usize) -> usize {
    let (row, col) = if row >= col { (row, col) } else { (col, row) };
    col * n_observed - col * (col.saturating_sub(1)) / 2 + (row - col)
}

fn add_symmetric_outer_vech(
    target: &mut nalgebra::DVectorViewMut<'_, f64>,
    left: &DMatrix<f64>,
    left_col: usize,
    left_support: &[usize],
    right: &DMatrix<f64>,
    right_col: usize,
    right_support: &[usize],
    n_observed: usize,
) {
    for &left_row in left_support {
        let left_value = left[(left_row, left_col)];
        for &right_row in right_support {
            let contribution = left_value * right[(right_row, right_col)];
            let index = vech_index(left_row, right_row, n_observed);
            target[index] += if left_row == right_row {
                2.0 * contribution
            } else {
                contribution
            };
        }
    }
}

fn add_self_outer_vech(
    target: &mut nalgebra::DVectorViewMut<'_, f64>,
    matrix: &DMatrix<f64>,
    col: usize,
    support: &[usize],
    n_observed: usize,
) {
    for (support_col, &observed_col) in support.iter().enumerate() {
        let col_value = matrix[(observed_col, col)];
        for &observed_row in &support[support_col..] {
            let index = vech_index(observed_row, observed_col, n_observed);
            target[index] += matrix[(observed_row, col)] * col_value;
        }
    }
}

fn solve_damped_normal_equations(
    normal_matrix: &DMatrix<f64>,
    rhs: &DVector<f64>,
) -> Option<DVector<f64>> {
    normal_matrix
        .clone()
        .cholesky()
        .map(|factor| factor.solve(rhs))
        .or_else(|| normal_matrix.clone().lu().solve(rhs))
}

fn invert_bread_matrix(bread_matrix: DMatrix<f64>) -> Option<DMatrix<f64>> {
    bread_matrix
        .clone()
        .cholesky()
        .map(|factor| factor.inverse())
        .or_else(|| bread_matrix.try_inverse())
}

fn ram_surfaces_from_rows(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    fixed_values: &[f64],
    free_values: &[f64],
    observed_indices: &[usize],
    free_row_offsets: &[i32],
    free_row_indices: &[i32],
    n_variables: usize,
) -> std::result::Result<(DMatrix<f64>, DMatrix<f64>), Error> {
    let (directed, covariance) = ram_matrices_from_rows(
        lhs_index,
        rhs_index,
        op_code,
        free_index,
        fixed_values,
        free_values,
        n_variables,
    );
    let (implied, implied_by_variable, observed_inverse) =
        implied_ram(&directed, &covariance, observed_indices)?;
    let delta = ram_delta_from_rows(
        lhs_index,
        rhs_index,
        op_code,
        free_row_offsets,
        free_row_indices,
        &implied_by_variable,
        &observed_inverse,
        free_values.len(),
    );

    Ok((implied, delta))
}

fn ram_implied_from_rows(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    fixed_values: &[f64],
    free_values: &[f64],
    observed_indices: &[usize],
    n_variables: usize,
) -> std::result::Result<DMatrix<f64>, Error> {
    let (directed, covariance) = ram_matrices_from_rows(
        lhs_index,
        rhs_index,
        op_code,
        free_index,
        fixed_values,
        free_values,
        n_variables,
    );
    let implied = implied_ram_observed(&directed, &covariance, observed_indices)?;

    Ok(implied)
}

fn fit_ram_dwls_with_plan(
    plan: &RamDwlsPlan,
    sample_cov_values: Vec<f64>,
    wls_values: Vec<f64>,
    initial_free_values: Vec<f64>,
    max_iter: i32,
    tol: f64,
    compute_se: bool,
) -> std::result::Result<List, Error> {
    if sample_cov_values.len() != plan.n_observed * plan.n_observed {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != plan.n_stats * plan.n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    if initial_free_values.len() != plan.n_free() {
        return Err(Error::Other("free_values has the wrong size".into()));
    }

    let sample_cov = from_col_major(&sample_cov_values, plan.n_observed, plan.n_observed);
    let wls_v = DwlsWeights::from_col_major(&wls_values, plan.n_stats);
    let observed = vech(&sample_cov);
    let mut params = DVector::from_vec(initial_free_values);
    for idx in 0..params.len() {
        params[idx] = params[idx]
            .max(plan.lower_bounds[idx])
            .min(plan.upper_bounds[idx]);
    }
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;
        let params_vec = params.as_slice();
        let (implied, delta) = ram_surfaces_from_rows(
            &plan.lhs_index,
            &plan.rhs_index,
            &plan.op_code,
            &plan.free_index,
            &plan.fixed_values,
            params_vec,
            &plan.observed_indices,
            &plan.free_row_offsets,
            &plan.free_row_indices,
            plan.n_variables,
        )?;
        let residual = &observed - vech(&implied);
        let weighted_residual = wls_v.apply_vector(&residual);
        let (gradient, hessian) = wls_v.normal_equations(&delta, &residual);
        let mut regularized = hessian.clone();

        for idx in 0..regularized.nrows() {
            regularized[(idx, idx)] += damping;
        }

        let Some(step) = solve_damped_normal_equations(&regularized, &gradient) else {
            damping *= 10.0;
            continue;
        };

        let old_objective = 0.5 * residual.dot(&weighted_residual);
        let mut next_params = &params + &step;

        for idx in 0..next_params.len() {
            next_params[idx] = next_params[idx]
                .max(plan.lower_bounds[idx])
                .min(plan.upper_bounds[idx]);
        }

        let next_implied = ram_implied_from_rows(
            &plan.lhs_index,
            &plan.rhs_index,
            &plan.op_code,
            &plan.free_index,
            &plan.fixed_values,
            next_params.as_slice(),
            &plan.observed_indices,
            plan.n_variables,
        )?;
        let next_residual = &observed - vech(&next_implied);
        let next_objective = wls_v.objective(&next_residual);

        if next_objective <= old_objective {
            params = next_params;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    let (implied, delta) = ram_surfaces_from_rows(
        &plan.lhs_index,
        &plan.rhs_index,
        &plan.op_code,
        &plan.free_index,
        &plan.fixed_values,
        params.as_slice(),
        &plan.observed_indices,
        &plan.free_row_offsets,
        &plan.free_row_indices,
        plan.n_variables,
    )?;
    let residual = &observed - vech(&implied);
    let naive_se = if compute_se {
        let bread = invert_bread_matrix(wls_v.weighted_crossprod(&delta))
            .unwrap_or_else(|| DMatrix::<f64>::zeros(params.len(), params.len()));
        bread.diagonal().map(|value| value.max(0.0).sqrt())
    } else {
        DVector::<f64>::zeros(params.len())
    };
    let objective = wls_v.objective(&residual);

    Ok(list!(
        estimates = params.as_slice().to_vec(),
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

/// Fit the covariance-only one-factor DWLS slice used by GenomicSEM's
/// `commonfactor()` model.
///
/// `sample_cov` and `wls_v` are flattened column-major R matrices.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param n Number of observed variables.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_one_factor_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    n: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if n <= 0 {
        return Err(Error::Other("n must be positive".into()));
    }

    let n = n as usize;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = DwlsWeights::from_col_major(&wls_values, n_stats);
    let observed = vech(&sample_cov);
    let (mut loadings, mut residuals) = initial_one_factor(&sample_cov);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;

        let implied = implied_one_factor(&loadings, &residuals);
        let residual = &observed - vech(&implied);
        let delta = delta_one_factor(&loadings);
        let weighted_residual = wls_v.apply_vector(&residual);
        let weighted_delta = wls_v.apply_matrix(&delta);
        let gradient = delta.transpose() * &weighted_residual;
        let hessian = delta.transpose() * &weighted_delta;
        let mut regularized = hessian.clone();

        for idx in 0..regularized.nrows() {
            regularized[(idx, idx)] += damping;
        }

        let Some(step) = solve_damped_normal_equations(&regularized, &gradient) else {
            damping *= 10.0;
            continue;
        };

        let old_objective = 0.5 * residual.dot(&weighted_residual);
        let mut next_loadings = loadings.clone();
        let mut next_residuals = residuals.clone();

        for idx in 0..n {
            next_loadings[idx] += step[idx];
            next_residuals[idx] = (next_residuals[idx] + step[n + idx]).max(1e-10);
        }

        let next_implied = implied_one_factor(&next_loadings, &next_residuals);
        let next_residual = &observed - vech(&next_implied);
        let next_objective = wls_v.objective(&next_residual);

        if next_objective <= old_objective {
            loadings = next_loadings;
            residuals = next_residuals;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    if loadings.sum() < 0.0 {
        loadings *= -1.0;
    }

    let implied = implied_one_factor(&loadings, &residuals);
    let residual = &observed - vech(&implied);
    let delta = delta_one_factor(&loadings);
    let weighted_delta = wls_v.apply_matrix(&delta);
    let bread = (delta.transpose() * &weighted_delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(2 * n, 2 * n));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = wls_v.objective(&residual);

    Ok(list!(
        loadings = loadings.as_slice().to_vec(),
        residuals = residuals.as_slice().to_vec(),
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

/// Fit a covariance-only observed-variable DWLS model where the free parameters
/// are selected entries of `vech(Sigma)`.
///
/// This covers the common-factor null model and the parameter-table refit used
/// to estimate its residual-covariance Q statistic.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param free_mask Integer mask over `vech(Sigma)`.
/// @param fixed_values Fixed values over `vech(Sigma)`.
/// @param n Number of observed variables.
/// @export
#[extendr]
fn fit_observed_covariance_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    free_mask: Integers,
    fixed_values: Doubles,
    n: i32,
) -> std::result::Result<List, Error> {
    if n <= 0 {
        return Err(Error::Other("n must be positive".into()));
    }

    let n = n as usize;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();
    let free_mask_values = free_mask
        .iter()
        .map(|value| value.0 != 0)
        .collect::<Vec<_>>();
    let fixed_values = fixed_values.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    if free_mask_values.len() != n_stats || fixed_values.len() != n_stats {
        return Err(Error::Other(
            "free_mask and fixed_values must match vech(sample_cov)".into(),
        ));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = DwlsWeights::from_col_major(&wls_values, n_stats);
    let observed = vech(&sample_cov);
    let fixed = DVector::from_vec(fixed_values);
    let free_indices = free_mask_values
        .iter()
        .enumerate()
        .filter_map(|(idx, is_free)| is_free.then_some(idx))
        .collect::<Vec<_>>();
    let mut delta = DMatrix::<f64>::zeros(n_stats, free_indices.len());

    for (col, row) in free_indices.iter().enumerate() {
        delta[(*row, col)] = 1.0;
    }

    let target = observed.clone() - fixed.clone();
    let weighted_delta = wls_v.apply_matrix(&delta);
    let weighted_target = wls_v.apply_vector(&target);
    let hessian = delta.transpose() * &weighted_delta;
    let rhs = delta.transpose() * &weighted_target;
    let Some(params) = hessian.clone().lu().solve(&rhs) else {
        return Err(Error::Other("observed covariance solve failed".into()));
    };

    let implied_vec = fixed + &delta * &params;
    let implied = from_vech(&implied_vec, n);
    let residual = observed - implied_vec;
    let bread = hessian
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(free_indices.len(), free_indices.len()));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = wls_v.objective(&residual);

    Ok(list!(
        estimates = params.as_slice().to_vec(),
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = true
    ))
}

/// Fit the marker-scaled common-factor GWAS DWLS slice used by
/// GenomicSEM's `commonfactorGWAS()` model.
///
/// Variable order is expected to be `SNP, trait_1, ..., trait_k`.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param k Number of trait indicators.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_commonfactor_gwas_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    k: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if k <= 0 {
        return Err(Error::Other("k must be positive".into()));
    }

    let k = k as usize;
    let n = k + 1;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = DwlsWeights::from_col_major(&wls_values, n_stats);
    let observed = vech(&sample_cov);
    let (mut loadings, mut gamma, mut residuals, mut psi, mut phi) =
        initial_commonfactor_gwas(&sample_cov);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;

        let implied = implied_commonfactor_gwas(&loadings, gamma, &residuals, psi, phi);
        let residual = &observed - vech(&implied);
        let delta = delta_commonfactor_gwas(&loadings, gamma, psi, phi);
        let weighted_residual = wls_v.apply_vector(&residual);
        let weighted_delta = wls_v.apply_matrix(&delta);
        let gradient = delta.transpose() * &weighted_residual;
        let hessian = delta.transpose() * &weighted_delta;
        let mut regularized = hessian.clone();

        for idx in 0..regularized.nrows() {
            regularized[(idx, idx)] += damping;
        }

        let Some(step) = regularized.lu().solve(&gradient) else {
            damping *= 10.0;
            continue;
        };

        let old_objective = 0.5 * residual.dot(&weighted_residual);
        let mut next_loadings = loadings.clone();
        let mut next_residuals = residuals.clone();

        for idx in 1..k {
            next_loadings[idx] += step[idx - 1];
        }

        let gamma_idx = k - 1;
        let theta_start = gamma_idx + 1;
        let next_gamma = gamma + step[gamma_idx];

        for idx in 0..k {
            next_residuals[idx] = (next_residuals[idx] + step[theta_start + idx]).max(1e-10);
        }

        let next_psi = (psi + step[theta_start + k]).max(1e-10);
        let next_phi = (phi + step[theta_start + k + 1]).max(1e-10);
        let next_implied = implied_commonfactor_gwas(
            &next_loadings,
            next_gamma,
            &next_residuals,
            next_psi,
            next_phi,
        );
        let next_residual = &observed - vech(&next_implied);
        let next_objective = wls_v.objective(&next_residual);

        if next_objective <= old_objective {
            loadings = next_loadings;
            gamma = next_gamma;
            residuals = next_residuals;
            psi = next_psi;
            phi = next_phi;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    let implied = implied_commonfactor_gwas(&loadings, gamma, &residuals, psi, phi);
    let residual = &observed - vech(&implied);
    let delta = delta_commonfactor_gwas(&loadings, gamma, psi, phi);
    let weighted_delta = wls_v.apply_matrix(&delta);
    let bread = (delta.transpose() * &weighted_delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(2 * k + 2, 2 * k + 2));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = wls_v.objective(&residual);

    Ok(list!(
        loadings = loadings.as_slice().to_vec(),
        gamma = gamma,
        residuals = residuals.as_slice().to_vec(),
        psi = psi,
        phi = phi,
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

/// Fit the common-factor GWAS Q-model refit where direct SNP effects and trait
/// residual variances are free and all first-stage quantities are fixed.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param loadings Fixed marker-scaled trait loadings.
/// @param gamma Fixed factor regression coefficient.
/// @param direct Initial direct SNP effects.
/// @param residuals Initial trait residual variances.
/// @param psi Fixed factor residual variance.
/// @param phi Fixed SNP variance.
/// @param k Number of trait indicators.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_commonfactor_gwas_q_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    loadings: Doubles,
    gamma: f64,
    direct: Doubles,
    residuals: Doubles,
    psi: f64,
    phi: f64,
    k: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if k <= 0 {
        return Err(Error::Other("k must be positive".into()));
    }

    let k = k as usize;
    let n = k + 1;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();
    let loadings = loadings.iter().map(|value| value.0).collect::<Vec<_>>();
    let direct = direct.iter().map(|value| value.0).collect::<Vec<_>>();
    let residuals = residuals.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    if loadings.len() != k || direct.len() != k || residuals.len() != k {
        return Err(Error::Other(
            "loadings, direct, and residuals must each have length k".into(),
        ));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = DwlsWeights::from_col_major(&wls_values, n_stats);
    let observed = vech(&sample_cov);
    let loadings = DVector::from_vec(loadings);
    let mut direct = DVector::from_vec(direct);
    let mut residuals = DVector::from_vec(residuals);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;

        let implied = implied_commonfactor_gwas_q(&loadings, gamma, &direct, &residuals, psi, phi);
        let residual = &observed - vech(&implied);
        let delta = delta_commonfactor_gwas_q(&loadings, gamma, &direct, phi);
        let weighted_residual = wls_v.apply_vector(&residual);
        let weighted_delta = wls_v.apply_matrix(&delta);
        let gradient = delta.transpose() * &weighted_residual;
        let hessian = delta.transpose() * &weighted_delta;
        let mut regularized = hessian.clone();

        for idx in 0..regularized.nrows() {
            regularized[(idx, idx)] += damping;
        }

        let Some(step) = regularized.lu().solve(&gradient) else {
            damping *= 10.0;
            continue;
        };

        let old_objective = 0.5 * residual.dot(&weighted_residual);
        let mut next_direct = direct.clone();
        let mut next_residuals = residuals.clone();

        for idx in 0..k {
            next_direct[idx] += step[idx];
            next_residuals[idx] = (next_residuals[idx] + step[k + idx]).max(1e-10);
        }

        let next_implied =
            implied_commonfactor_gwas_q(&loadings, gamma, &next_direct, &next_residuals, psi, phi);
        let next_residual = &observed - vech(&next_implied);
        let next_objective = wls_v.objective(&next_residual);

        if next_objective <= old_objective {
            direct = next_direct;
            residuals = next_residuals;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    let implied = implied_commonfactor_gwas_q(&loadings, gamma, &direct, &residuals, psi, phi);
    let residual = &observed - vech(&implied);
    let delta = delta_commonfactor_gwas_q(&loadings, gamma, &direct, phi);
    let weighted_delta = wls_v.apply_matrix(&delta);
    let bread = (delta.transpose() * &weighted_delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(2 * k, 2 * k));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = wls_v.objective(&residual);

    Ok(list!(
        direct = direct.as_slice().to_vec(),
        residuals = residuals.as_slice().to_vec(),
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

/// Fit the fixed-measurement one-factor GWAS DWLS slice used by the default
/// `userGWAS(..., fix_measurement = TRUE)` path.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param loadings Fixed marker-scaled trait loadings.
/// @param residuals Initial trait residual variances.
/// @param psi Initial factor residual variance.
/// @param gamma Initial factor regression coefficient.
/// @param phi Initial SNP variance.
/// @param k Number of trait indicators.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_user_gwas_fixed_measurement_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    loadings: Doubles,
    residuals: Doubles,
    psi: f64,
    gamma: f64,
    phi: f64,
    k: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if k <= 0 {
        return Err(Error::Other("k must be positive".into()));
    }

    let k = k as usize;
    let n = k + 1;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();
    let loadings = loadings.iter().map(|value| value.0).collect::<Vec<_>>();
    let residuals = residuals.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    if loadings.len() != k || residuals.len() != k {
        return Err(Error::Other(
            "loadings and residuals must each have length k".into(),
        ));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = DwlsWeights::from_col_major(&wls_values, n_stats);
    let observed = vech(&sample_cov);
    let loadings = DVector::from_vec(loadings);
    let mut residuals = DVector::from_vec(residuals);
    let mut psi = psi.max(1e-10);
    let mut gamma = gamma;
    let mut phi = phi.max(1e-10);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;

        let implied = implied_commonfactor_gwas(&loadings, gamma, &residuals, psi, phi);
        let residual = &observed - vech(&implied);
        let delta = delta_user_gwas_fixed_measurement(&loadings, gamma, phi);
        let weighted_residual = wls_v.apply_vector(&residual);
        let weighted_delta = wls_v.apply_matrix(&delta);
        let gradient = delta.transpose() * &weighted_residual;
        let hessian = delta.transpose() * &weighted_delta;
        let mut regularized = hessian.clone();

        for idx in 0..regularized.nrows() {
            regularized[(idx, idx)] += damping;
        }

        let Some(step) = regularized.lu().solve(&gradient) else {
            damping *= 10.0;
            continue;
        };

        let old_objective = 0.5 * residual.dot(&weighted_residual);
        let mut next_residuals = residuals.clone();

        for idx in 0..k {
            next_residuals[idx] = (next_residuals[idx] + step[idx]).max(1e-10);
        }

        let next_psi = (psi + step[k]).max(1e-10);
        let next_gamma = gamma + step[k + 1];
        let next_phi = (phi + step[k + 2]).max(1e-10);
        let next_implied =
            implied_commonfactor_gwas(&loadings, next_gamma, &next_residuals, next_psi, next_phi);
        let next_residual = &observed - vech(&next_implied);
        let next_objective = wls_v.objective(&next_residual);

        if next_objective <= old_objective {
            residuals = next_residuals;
            psi = next_psi;
            gamma = next_gamma;
            phi = next_phi;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    let implied = implied_commonfactor_gwas(&loadings, gamma, &residuals, psi, phi);
    let residual = &observed - vech(&implied);
    let delta = delta_user_gwas_fixed_measurement(&loadings, gamma, phi);
    let weighted_delta = wls_v.apply_matrix(&delta);
    let bread = (delta.transpose() * &weighted_delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(k + 3, k + 3));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = wls_v.objective(&residual);

    Ok(list!(
        residuals = residuals.as_slice().to_vec(),
        psi = psi,
        gamma = gamma,
        phi = phi,
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

/// Evaluate generic RAM-model implied covariance and Jacobian surfaces.
///
/// Row indices are 1-based to match R. Operator codes are:
/// `1 = =~`, `2 = ~`, `3 = ~~`.
/// @param lhs_index Parameter-table lhs row indices.
/// @param rhs_index Parameter-table rhs row indices.
/// @param op_code Parameter-table operator codes.
/// @param free_index Free-parameter indices, with `0` for fixed rows.
/// @param fixed_values Fixed row values, ignored for free rows.
/// @param free_values Current free-parameter values.
/// @param observed_index Observed-variable indices in the full RAM system.
/// @param free_row_offsets Offsets into `free_row_indices` for each free parameter.
/// @param free_row_indices Flattened 1-based row indices grouped by free parameter.
/// @param n_variables Number of variables in the full RAM system.
/// @export
#[extendr]
fn evaluate_ram_surfaces(
    lhs_index: Integers,
    rhs_index: Integers,
    op_code: Integers,
    free_index: Integers,
    fixed_values: Doubles,
    free_values: Doubles,
    observed_index: Integers,
    free_row_offsets: Integers,
    free_row_indices: Integers,
    n_variables: i32,
) -> std::result::Result<List, Error> {
    if n_variables <= 0 {
        return Err(Error::Other("n_variables must be positive".into()));
    }

    let n_variables = n_variables as usize;
    let lhs_index = lhs_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let rhs_index = rhs_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let op_code = op_code.iter().map(|value| value.0).collect::<Vec<_>>();
    let free_index = free_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let fixed_values = fixed_values.iter().map(|value| value.0).collect::<Vec<_>>();
    let free_values = free_values.iter().map(|value| value.0).collect::<Vec<_>>();
    let observed_index = observed_index
        .iter()
        .map(|value| value.0)
        .collect::<Vec<_>>();
    let free_row_offsets = free_row_offsets
        .iter()
        .map(|value| value.0)
        .collect::<Vec<_>>();
    let free_row_indices = free_row_indices
        .iter()
        .map(|value| value.0)
        .collect::<Vec<_>>();

    let observed_indices = validate_ram_indices(
        &lhs_index,
        &rhs_index,
        &op_code,
        &free_index,
        &fixed_values,
        &observed_index,
        n_variables,
        free_values.len(),
    )?;
    validate_free_row_groups(
        &free_row_offsets,
        &free_row_indices,
        &free_index,
        free_values.len(),
    )?;
    let (implied, delta) = ram_surfaces_from_rows(
        &lhs_index,
        &rhs_index,
        &op_code,
        &free_index,
        &fixed_values,
        &free_values,
        &observed_indices,
        &free_row_offsets,
        &free_row_indices,
        n_variables,
    )?;

    Ok(list!(
        implied = to_col_major(&implied),
        delta = to_col_major(&delta)
    ))
}

/// Compile immutable generic RAM-model structure into a reusable native plan.
/// @param lhs_index Parameter-table lhs row indices.
/// @param rhs_index Parameter-table rhs row indices.
/// @param op_code Parameter-table operator codes.
/// @param free_index Free-parameter indices, with `0` for fixed rows.
/// @param fixed_values Fixed row values, ignored for free rows.
/// @param observed_index Observed-variable indices in the full RAM system.
/// @param free_row_offsets Offsets into `free_row_indices` for each free parameter.
/// @param free_row_indices Flattened 1-based row indices grouped by free parameter.
/// @param lower_bounds Lower bound for each free parameter.
/// @param upper_bounds Upper bound for each free parameter.
/// @param n_variables Number of variables in the full RAM system.
/// @export
#[extendr]
fn compile_ram_dwls_plan(
    lhs_index: Integers,
    rhs_index: Integers,
    op_code: Integers,
    free_index: Integers,
    fixed_values: Doubles,
    observed_index: Integers,
    free_row_offsets: Integers,
    free_row_indices: Integers,
    lower_bounds: Doubles,
    upper_bounds: Doubles,
    n_variables: i32,
) -> std::result::Result<Robj, Error> {
    if n_variables <= 0 {
        return Err(Error::Other("n_variables must be positive".into()));
    }

    let plan = RamDwlsPlan::new(
        lhs_index.iter().map(|value| value.0).collect::<Vec<_>>(),
        rhs_index.iter().map(|value| value.0).collect::<Vec<_>>(),
        op_code.iter().map(|value| value.0).collect::<Vec<_>>(),
        free_index.iter().map(|value| value.0).collect::<Vec<_>>(),
        fixed_values.iter().map(|value| value.0).collect::<Vec<_>>(),
        observed_index
            .iter()
            .map(|value| value.0)
            .collect::<Vec<_>>(),
        free_row_offsets
            .iter()
            .map(|value| value.0)
            .collect::<Vec<_>>(),
        free_row_indices
            .iter()
            .map(|value| value.0)
            .collect::<Vec<_>>(),
        lower_bounds.iter().map(|value| value.0).collect::<Vec<_>>(),
        upper_bounds.iter().map(|value| value.0).collect::<Vec<_>>(),
        n_variables as usize,
    )?;

    Ok(ExternalPtr::new(plan).into())
}

/// Check whether a generic RAM native plan still points at live Rust state.
///
/// External pointers are intentionally not serializable across all R worker
/// strategies, so the R layer can rebuild them lazily when needed.
/// @param plan Native plan external pointer.
/// @export
#[extendr]
fn ram_dwls_plan_is_valid(plan: Robj) -> bool {
    ExternalPtr::<RamDwlsPlan>::try_from(plan)
        .and_then(|plan| plan.try_addr().map(|_| ()))
        .is_ok()
}

/// Fit a generic RAM-model DWLS slice from a reusable native plan.
///
/// @param plan Native plan external pointer from `compile_ram_dwls_plan()`.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param free_values Initial free-parameter values.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @param compute_se Whether to compute naive standard errors from the final bread matrix.
/// @export
#[extendr]
fn fit_ram_dwls_plan(
    plan: Robj,
    sample_cov: Doubles,
    wls_v: Doubles,
    free_values: Doubles,
    max_iter: i32,
    tol: f64,
    compute_se: bool,
) -> std::result::Result<List, Error> {
    let plan: ExternalPtr<RamDwlsPlan> = plan.try_into()?;
    let plan = plan.try_addr()?;

    fit_ram_dwls_with_plan(
        plan,
        sample_cov.iter().map(|value| value.0).collect::<Vec<_>>(),
        wls_v.iter().map(|value| value.0).collect::<Vec<_>>(),
        free_values.iter().map(|value| value.0).collect::<Vec<_>>(),
        max_iter,
        tol,
        compute_se,
    )
}

/// Fit a generic RAM-model DWLS slice from compiled row data.
///
/// This compatibility entry point rebuilds a native plan for one-off callers.
/// Repeated fits should prefer `compile_ram_dwls_plan()` plus
/// `fit_ram_dwls_plan()`.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param lhs_index Parameter-table lhs row indices.
/// @param rhs_index Parameter-table rhs row indices.
/// @param op_code Parameter-table operator codes.
/// @param free_index Free-parameter indices, with `0` for fixed rows.
/// @param fixed_values Fixed row values, ignored for free rows.
/// @param free_values Initial free-parameter values.
/// @param observed_index Observed-variable indices in the full RAM system.
/// @param free_row_offsets Offsets into `free_row_indices` for each free parameter.
/// @param free_row_indices Flattened 1-based row indices grouped by free parameter.
/// @param lower_bounds Lower bound for each free parameter.
/// @param upper_bounds Upper bound for each free parameter.
/// @param n_variables Number of variables in the full RAM system.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @param compute_se Whether to compute naive standard errors from the final bread matrix.
/// @export
#[extendr]
fn fit_ram_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    lhs_index: Integers,
    rhs_index: Integers,
    op_code: Integers,
    free_index: Integers,
    fixed_values: Doubles,
    free_values: Doubles,
    observed_index: Integers,
    free_row_offsets: Integers,
    free_row_indices: Integers,
    lower_bounds: Doubles,
    upper_bounds: Doubles,
    n_variables: i32,
    max_iter: i32,
    tol: f64,
    compute_se: bool,
) -> std::result::Result<List, Error> {
    if n_variables <= 0 {
        return Err(Error::Other("n_variables must be positive".into()));
    }

    let plan = RamDwlsPlan::new(
        lhs_index.iter().map(|value| value.0).collect::<Vec<_>>(),
        rhs_index.iter().map(|value| value.0).collect::<Vec<_>>(),
        op_code.iter().map(|value| value.0).collect::<Vec<_>>(),
        free_index.iter().map(|value| value.0).collect::<Vec<_>>(),
        fixed_values.iter().map(|value| value.0).collect::<Vec<_>>(),
        observed_index
            .iter()
            .map(|value| value.0)
            .collect::<Vec<_>>(),
        free_row_offsets
            .iter()
            .map(|value| value.0)
            .collect::<Vec<_>>(),
        free_row_indices
            .iter()
            .map(|value| value.0)
            .collect::<Vec<_>>(),
        lower_bounds.iter().map(|value| value.0).collect::<Vec<_>>(),
        upper_bounds.iter().map(|value| value.0).collect::<Vec<_>>(),
        n_variables as usize,
    )?;

    fit_ram_dwls_with_plan(
        &plan,
        sample_cov.iter().map(|value| value.0).collect::<Vec<_>>(),
        wls_v.iter().map(|value| value.0).collect::<Vec<_>>(),
        free_values.iter().map(|value| value.0).collect::<Vec<_>>(),
        max_iter,
        tol,
        compute_se,
    )
}

extendr_module! {
    mod lavaanrust;
    fn fit_one_factor_dwls;
    fn fit_observed_covariance_dwls;
    fn fit_commonfactor_gwas_dwls;
    fn fit_commonfactor_gwas_q_dwls;
    fn fit_user_gwas_fixed_measurement_dwls;
    fn evaluate_ram_surfaces;
    fn compile_ram_dwls_plan;
    fn ram_dwls_plan_is_valid;
    fn fit_ram_dwls_plan;
    fn fit_ram_dwls;
}
