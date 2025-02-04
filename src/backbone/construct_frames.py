

"""
def get_vector_info_single(self, basis_vectors: chex.Array) -> Tuple[chex.Array, chex.Array, chex.Array]:
    chex.assert_rank(basis_vectors, 2)
    n_vectors, dim = basis_vectors.shape
    assert dim == 3
    assert n_vectors == 2 or n_vectors == 3
    basis_vectors = basis_vectors
    vec1 = basis_vectors[0]
    vec2 = basis_vectors[1]
    arccos_in = jnp.dot(vec1, vec2) / safe_norm(vec1, axis=-1) / safe_norm(vec2, axis=-1)
    theta = jnp.arccos(arccos_in)
    log_barrier_in = 1 - jnp.abs(arccos_in)

    # Don't penalize when very close (or on top of each other) to keep grads stable.
    log_barrier_in = jnp.where(log_barrier_in < 1e-4, jnp.ones_like(log_barrier_in), log_barrier_in)
    aux_loss = - jnp.log(log_barrier_in)
    return theta, aux_loss, log_barrier_in
    """