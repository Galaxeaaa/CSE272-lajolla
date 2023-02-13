#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
        dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Vector3 h = (dir_in + dir_out) / length(dir_in + dir_out);
    Vector3 h_local = to_local(frame, h);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real h_dot_out = dot(h, dir_out);
    Real n_dot_in = dot(frame.n, dir_in);

    Spectrum F = schlick_fresnel(base_color, abs(h_dot_out));
    Real D = GTR2(h_local, roughness, anisotropic);
    Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
             smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);

    Spectrum f_metal = (F * D * G) / (4 * abs(n_dot_in));

    return f_metal;
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
        dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Vector3 h = (dir_in + dir_out) / length(dir_in + dir_out);
    Vector3 h_local = to_local(frame, h);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real h_dot_out = dot(h, dir_out);
    Real n_dot_in = dot(frame.n, dir_in);

    Real D = GTR2(h_local, roughness, anisotropic);
    Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic);

    return (D * G) / (4 * abs(n_dot_in));
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 local_micro_normal = sample_visible_normals(to_local(frame, dir_in), alpha_x, alpha_y, rnd_param_uv);

    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
