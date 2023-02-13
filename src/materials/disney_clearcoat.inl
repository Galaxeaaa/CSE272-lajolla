#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;

    Real h_dot_out = dot(h, dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);

    auto R_0 = [](Real eta) {
        return ((eta - 1) * (eta - 1)) / ((eta + 1) * (eta + 1));
    };

    Real F = R_0(1.5) + (1 - R_0(1.5)) * pow((1 - abs(h_dot_out)), 5);

    Real D = (alpha * alpha - 1) / (c_PI * log(alpha * alpha) * (1 + (alpha * alpha - 1) * h_local.z * h_local.z));

    auto Lambda = [&](Vector3 dir) {
        Vector3 dir_local = to_local(frame, dir);
        return (sqrt(1 + (((dir_local.x * 0.25) * (dir_local.x * 0.25) + (dir_local.y * 0.25) * (dir_local.y * 0.25)) /
                          (dir_local.z * dir_local.z))) -
                1) /
               2;
    };

    auto G_ = [&](Vector3 dir) {
        return 1 / (1 + Lambda(dir));
    };

    Real G = G_(dir_in) * G_(dir_out);

    Real clearcoat = (F * D * G) / (4 * abs(n_dot_in));
    Spectrum f_clearcoat = Spectrum(clearcoat, clearcoat, clearcoat);

    return f_clearcoat;
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;

    Real h_dot_out = dot(h, dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_h = dot(frame.n, h);

    Real D = (alpha * alpha - 1) / (c_PI * log(alpha * alpha) * (1 + (alpha * alpha - 1) * h_local.z * h_local.z));

    return (D * n_h) / (4 * abs(h_dot_out));
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;

    Real cos_h_elevation = sqrt((1 - pow(alpha * alpha, 1 - rnd_param_uv.x)) / (1 - alpha * alpha));
    Real sin_h_elevation = sqrt(1 - cos_h_elevation * cos_h_elevation);
    Real h_azimuth = c_TWOPI * rnd_param_uv.y;

    Vector3 h_local = Vector3(sin_h_elevation * cos(h_azimuth), sin_h_elevation * sin(h_azimuth), cos_h_elevation);
    Vector3 h = to_world(frame, h_local);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, Real(1) /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
