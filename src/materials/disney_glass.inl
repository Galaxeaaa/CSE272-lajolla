#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) * dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
    }
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }
    Vector3 h_local = to_local(frame, h);
    Real h_dot_in = dot(h, dir_in);
    Real h_dot_out = dot(h, dir_out);
    Real n_dot_in = dot(frame.n, dir_in);

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));

    Real F = fresnel_dielectric(h_dot_in, eta);
    // // Schlick
    // Real F_0 = pow((1 - eta) / (1 + eta), 2);
    // Real cos_theta_t_sq = 1 - (1 - pow(cos(n_dot_in), 2)) / (eta * eta);
    // Real F;
    // if (cos_theta_t_sq > 0) {
    //     F = F_0 + (1 - F_0) * pow(1 - ndot, 5);
    // }
    // else {
    //     F = 1;
    // }

    Real D = GTR2(h_local, roughness, anisotropic);
    Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
             smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);

    if (reflect) {
        return base_color * (F * D * G) / (4 * fabs(n_dot_in));
    } else {
        Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
        return sqrt(base_color) * eta_factor * eta * eta * (1 - F) * D * G * fabs(h_dot_out * h_dot_in) / (fabs(n_dot_in) * (h_dot_in + eta * h_dot_out) * (h_dot_in + eta * h_dot_out));
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                       dot(vertex.geometric_normal, dir_out) >
                   0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }
    Vector3 h_local = to_local(frame, h);

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(h, dir_in);
    Real h_dot_out = dot(h, dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
	anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));

    Real F = fresnel_dielectric(h_dot_in, eta);
    Real D = GTR2(h_local, roughness, anisotropic);
    Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic);

    if (reflect) {
        return (F * D * G) / (4 * fabs(n_dot_in));
    } else {
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F) * D * G * fabs(dh_dout * h_dot_in / n_dot_in);
    }
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anistropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aspect = sqrt(1 - 0.9 * anistropic);
    Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 local_micro_normal = sample_visible_normals(to_local(frame, dir_in), alpha_x, alpha_y, rnd_param_uv);

    // Transform the micro normal to world space
    Vector3 h = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(h, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half vector if needed
        if (h_dot_in < 0) {
            h = -h;
        }
        Real h_dot_out = sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * h;
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
