# Single monochromatic absorption-only homogeneous volume

1. I see a circle filled with uniform light blue color. Because the volume doesn't absorb any light, so we can see the sphere as if there's nothing between.

2. The solution to the ODE
    $$\frac{d}{dt}L(p(t), \omega) = -\sigma_aL(p(t), \omega) + L_e(p(t), \omega)$$
    is
    $$L(p(t), \omega) = \int_0^t{T(t')L_e(p(t'), \omega)dt'}$$
    , where $T(t') = \exp\left(-\int_0^{t'}\sigma_adt''\right) = \exp(-\sigma_at)$. To estimate the integration, we can sample the transmittance in the same way with the following questions. If $t < t_{hit}$, we add $\frac{TL_e}{pdf}$ to the result. If $t \geq t_{hit}$, we apply the same pdf as the following questions and add the contribution of the emission of the surface to the result.
    ```python
    def L(screen_pos, rng):
        ...
        isect = intersect(scene, camera_ray)
        u = next(rng)
        t = -log(1 - u) / sigma_a
        if t < isect.t_hit:
            trans_pdf = exp(-sigma_a * t) * sigma_a
            trans = exp(-sigma_a * t)
        else:
            trans_pdf = exp(-sigma_a * t)
            trans = exp(-sigma_a * t)
        return trans * Le(t) / trans_pdf
        ...

    ```

# Single monochromatic homogeneous volume with absorption and single-scattering, no surface lighting

1. Since $p(t) \propto \exp(-\sigma_tt)$, we let $p(t) = C\exp(-\sigma_tt)$ and integrate it:
    $$P(t < t_{hit}) = \int_0^{t_{hit}}{C\exp(-\sigma_tt)dt} = C(-\frac{\exp(-\sigma_tt_{hit})}{\sigma_t}+\frac{1}{\sigma_t}) < \frac{C}{\sigma_t}$$
    When $t_{hit} \rightarrow \infty$, $P(t < t_{hit}) \rightarrow \frac{C}{\sigma_t}$. We want to sample $t$ in $[0, t_{hit}]$, so the natural choice of $C$ is $\sigma_t$, hence $p(t) = \sigma_t\exp(-\sigma_tt)$.

2. When $t \geq t_{hit}$, the ray hit the surface and we account for the surface emission with the transmittance pdf being $\exp(-\sigma_tt_{hit})$. It is done this way because if $t \geq t_{hit}$, the sampled point would go through the surface, which is not reasonable in the context of this question, hence we consider the ray hit the surface when $t \geq t_{hit}$ to make full use of the sampled distance.

3. The larger the $\sigma_a$, the darker the whole image, because the medium would absorb more radiance. The larger the $\sigma_s$, the more the lighting would spread out, because the medium is better at scattering light.

4. $g$ is the anisotropic factor. Large $g$ means the phase function value is large at $0\degree$ direction and small at $180\degree$ direction. So the larger the $g$, the more possible light would scatter at the same direction of its propagation. In this test scene, when we increase $g$, the radiance will concentrate more near the region of the light source.

# Multiple monochromatic homogeneous volumes with absorption and multiple-scattering using only phase function sampling, no surface lighting

1. For the volume filling the space, if we increase/decrease $\sigma_s$, the blue light will scatter in a larger/smaller scale; if we increase/decrease $\sigma_a$, the whole scene will appear darker/brighter. For the sphere volume, if we increase/decrease $\sigma_s$, the sphere light will appear less/more transparent throught the volume; if we increase/decrease $\sigma_a$, the volume will appear darker/brighter. If we increase/decrease `max_depth`, the image will be brighter/darker.

    If $\sigma_s$ and $\sigma_a$ are large, we can set a lower `max_depth`, since the radiance will, in this case, have a larger probability to be reduced to zero before reaching the `max_depth`.

2. When we increase $g$, the radiance will concentrate more near the region of the light source, out of the same reason with the last section. Additionally, the contrast of the sphere volume and the background volume becomes smaller.

3. The phase function value should be the same when the angle of input direction and output direction are the same. Given the input direction, the integral of the phase function over the sphere domain should be $1$ to be conservative. I think the phase function can be a guassian-like function, with two parameters controlling the mean and variance.