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

5. ***I did the bonus question.*** I implemented equiangular sampling with multi-importance sampling. The output image is shown in `vol2_bonus.exr`.

# Multiple monochromatic homogeneous volumes with absorption and multiple-scattering using only phase function sampling, no surface lighting

1. For the volume filling the space, if we increase/decrease $\sigma_s$, the blue light will scatter in a larger/smaller scale; if we increase/decrease $\sigma_a$, the whole scene will appear darker/brighter. For the sphere volume, if we increase/decrease $\sigma_s$, the sphere light will appear less/more transparent throught the volume; if we increase/decrease $\sigma_a$, the volume will appear darker/brighter. If we increase/decrease `max_depth`, the image will be brighter/darker.

    If $\sigma_s$ and $\sigma_a$ are large, we can set a lower `max_depth`, since the radiance will, in this case, have a larger probability to be reduced to zero before reaching the `max_depth`.

2. When we increase $g$, the radiance will concentrate more near the region of the light source, out of the same reason with the last section. Additionally, the contrast of the sphere volume and the background volume becomes smaller.

3. The phase function value should be the same when the angle of input direction and output direction are the same. Given the input direction, the integral of the phase function over the sphere domain should be $1$ to be conservative. I think the phase function can be a guassian-like function, with two parameters controlling the mean and variance.

# Multiple monochromatic homogeneous volumes with absorption and multiple-scattering with both phase function sampling and next event estimation, no surface lighting

1. When light source is large and it has a higher probability for a sampled light point to be reachable, NEE will be more efficient than phase function sampling. In the test scenes, the first scene will be more efficient using NEE, since it has a much smaller light source.

2. The radiance on the surface of an object composed of dense volume looks more evenly distributed compared to a lambertian surface. This is because the phase function of the dense volume is isotropic, while the lambertian surface only reflect the light to a specific direction.

3. Volume representation of scenes is a more general case compared to mesh representation and many other representations, and this is why Jim Kajiya stated like that in 1991. The reason why it hasn't happened yet is that volume rendering is too computationally expensive.

# Multiple monochromatic homogeneous volumes with absorption and multiple-scattering with both phase function sampling and next event estimation, with surface lighting

1. When we increase the index of refraction of the dielectric interface, the radiance from the light source on the other side passing through the dielectric sphere becomes more concentrated. This is because when we increase the index, the angle change of a ray passing through the interface becomes bigger, and the convex lens effect becomes stronger.

2. With blue interior medium and transparent glass, the handle of the lid appears darker. With blue glass and empty medium, the handle of the lid appears more shiny. Maybe this is because the teapot model has a filled handle of the lid. But anyway, this is a difference. Also, the specularity from the back is more obvious on the second case, since there's no interior medium scattering the specular light.

# Multiple chromatic homogeneous volumes with absorption and multiple-scattering with both phase function sampling and next event estimation, with surface lighting

1. For a heterogeneous volume with only a few places having large density, the null scattering is inefficient, because there will be many null collisions and many rejections, resulting in expensive construction of free-path samples.
    
    In this case, we can set a non-bounding majority $\sigma_m$ covering most of the volume. Previously we set the probability of sampling null particles as $\frac{\sigma_n}{\sigma_m}$ and that of sampling real particles as $\frac{\sigma_t}{\sigma_m}$, where $\sigma_n + \sigma_t = \sigma_m$. Now we don't use the majority to define the pdf but instead set the probability of sampling null particles as $\frac{|\sigma_n|}{\sigma_t + |\sigma_n|}$ and that of sampling real particles as $\frac{\sigma_a + \sigma_s}{\sigma_t + |\sigma_n|}$. This makes the algorithm more robust and be able to handle the previously inefficient case.

2. With emissive volume, the radiative transfer equation is as follows:

    $$\frac{d}{dt}L(\mathbf{p}(t), \omega) = -\sigma_mL(\mathbf{p}(t), \omega) + L_e(\mathbf{p}(t), \omega) + \sigma_n(\mathbf{p}(t))L(\mathbf{p}(t), \omega) + \sigma_s(\mathbf{p}(t))\int_{S^2}{\rho(\omega, \omega')L(\mathbf{p}(t), \omega')d\omega'}$$

    The solution of this ODE is

    $$\begin{aligned} L(p(t), \omega) = &\int_0^{t_{hit}}{T_m(p(0), p(t'))\left(L_e(\mathbf{p}(t'), \omega) + \sigma_n(\mathbf{p}(t'))L(\mathbf{p}(t'), \omega) + \sigma_s(\mathbf{p}(t'))\int_{S^2}{\rho(\omega, \omega')L(\mathbf{p}(t'), \omega')d\omega'}\right)dt'} \\
    &+ T_m(p(0), p(t_{hit}))L_{surface}(p(t_{hit})) \end{aligned}$$

    After we sample the homogenized transmittance $T_m$ and get a distance $t'$, if we hit a *real* particle, we not only compute the $S^2$ integral via importance sampling but also add the particle emission to the contribution in this iteration.

3. Unbiased solution makes the rendering physically accurate. Generally biased volume rendering solution is not a good idea, but for real-time rendering of relatively simple volumetric scenes, I think it is acceptable to use biased but fast volumetric rendering to just achieve aesthetic images. For example, we can assume uniform volume, or use lower resolution, or use smaller sampling rate to reduce computation complexity.