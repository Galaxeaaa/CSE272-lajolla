#include "../microfacet.h"

// Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
//     bool reflect = dot(vertex.geometric_normal, dir_in) *
//                        dot(vertex.geometric_normal, dir_out) >
//                    0;
//     // Flip the shading frame if it is inconsistent with the geometry normal
//     Frame frame = vertex.shading_frame;
//     if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
//         frame = -frame;
//     }
//     Real t = dot(vertex.geometric_normal, dir_in);
//     Real t2 = dot(vertex.geometric_normal, dir_out);

//     // Homework 1: implement this!
//     // all parameters required
//     Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
//     Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
//     // clamp
//     specular_transmission = std::clamp(specular_transmission, Real(0.01), Real(1));
//     metallic = std::clamp(metallic, Real(0.01), Real(1));
//     subsurface = std::clamp(subsurface, Real(0.01), Real(1));
//     specular = std::clamp(specular, Real(0.01), Real(1));
//     roughness = std::clamp(roughness, Real(0.01), Real(1));
//     specular_tint = std::clamp(specular_tint, Real(0.01), Real(1));
//     anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
//     sheen = std::clamp(sheen, Real(0.01), Real(1));
//     sheen_tint = std::clamp(sheen_tint, Real(0.01), Real(1));
//     clearcoat = std::clamp(clearcoat, Real(0.01), Real(1));
//     clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));
//     // tmp vector required
//     Vector3 h, h_local, local_dir_in, local_dir_out;

//     if (reflect) {
//         h = normalize(dir_in + dir_out);
//     } else {
//         h = normalize(dir_in + dir_out * eta);
//     }
//     if (dot(h, frame.n) < 0) {
//         h = -h;
//     }
//     h_local = to_local(frame, h);
//     local_dir_in = to_local(frame, dir_in);
//     local_dir_out = to_local(frame, dir_out);

//     Real h_dot_in = dot(h, dir_in);
//     Real h_dot_out = dot(h, dir_out);
//     Real n_dot_in = dot(frame.n, dir_in);
//     Real n_dot_out = dot(frame.n, dir_out);

//     Spectrum f_diffuse, f_metal, f_glass, f_clearcoat, f_sheen;
//     Real diffuseWeight, metalWeight, glassWeight, clearcoatWeight;
//     f_diffuse = f_metal = f_glass = f_clearcoat = f_sheen = make_zero_spectrum();
//     diffuseWeight = (1 - metallic) * (1 - specular_transmission);
//     metalWeight = (1 - specular_transmission * (1 - metallic));
//     glassWeight = (1 - metallic) * specular_transmission;
//     clearcoatWeight = 0.25 * clearcoat;

//     auto R_0 = [](Real eta) {
//         return ((eta - 1) * (eta - 1)) / ((eta + 1) * (eta + 1));
//     };

//     if (dot(vertex.geometric_normal, dir_in) > 0) {
//         if (dot(vertex.geometric_normal, dir_out) > 0) {
//             {
//                 // diffuse
//                 Real FD90 = 0.5 + 2 * roughness * h_dot_out * h_dot_out;
//                 Real FD_in = 1 + (FD90 - 1) * pow(1 - abs(n_dot_in), 5);
//                 Real FD_out = 1 + (FD90 - 1) * pow(1 - abs(n_dot_in), 5);
//                 Spectrum f_base = (base_color / c_PI) * FD_in * FD_out * n_dot_out;

//                 Real FSS90 = roughness * h_dot_out * h_dot_out;
//                 Real FSS_in = 1 + (FSS90 - 1) * pow(1 - abs(n_dot_in), 5);
//                 Real FSS_out = 1 + (FSS90 - 1) * pow(1 - abs(n_dot_out), 5);
//                 Spectrum f_subsurface = (1.25 * base_color / c_PI) * (FSS_in * FSS_out * (1 / (abs(n_dot_in) + abs(n_dot_out)) - 0.5) + 0.5) * abs(n_dot_out);

//                 f_diffuse = (1 - subsurface) * f_base + subsurface * f_subsurface;
//             }
//         }
//         {
//             // sheen
//             Spectrum C_tint;
//             if (luminance(base_color) > 0) {
//                 C_tint = base_color / luminance(base_color);
//             } else {
//                 C_tint = Spectrum(1, 1, 1);
//             }

//             Spectrum C_sheen = (1 - sheen_tint) + sheen_tint * C_tint;

//             f_sheen = C_sheen * pow(1 - abs(h_dot_out), 5) * abs(n_dot_out);

//             // metal
//             Spectrum K_s = (1 - specular_tint) * Spectrum(1, 1, 1) + specular_tint * C_tint;
//             Spectrum C_0 = specular * R_0(eta) * (1 - metallic) * K_s + metallic * base_color;
//             Spectrum F = schlick_fresnel(C_0, h_dot_out);

//             Real D = GTR2(h_local, roughness, anisotropic);

//             Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
//                      smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);

//             f_metal = (F * D * G) / (4 * abs(n_dot_in));
//         }
//         {
//             // clearcoat
//             Real alpha = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
//             Real F = R_0(1.5) + (1 - R_0(1.5)) * pow((1 - abs(h_dot_out)), 5);

//             Real D = (alpha * alpha - 1) / (c_PI * log(alpha * alpha) * (1 + (alpha * alpha - 1) * h_local.z * h_local.z));

//             auto Lambda = [&](Vector3 dir) {
//                 Vector3 dir_local = to_local(frame, dir);
//                 return (sqrt(1 + (((dir_local.x * 0.25) * (dir_local.x * 0.25) + (dir_local.y * 0.25) * (dir_local.y * 0.25)) /
//                                   (dir_local.z * dir_local.z))) -
//                         1) /
//                        2;
//             };

//             auto G_ = [&](Vector3 dir) {
//                 return 1 / (1 + Lambda(dir));
//             };

//             Real G = G_(dir_in) * G_(dir_out);

//             Real clearcoat = (F * D * G) / (4 * abs(n_dot_in));
//             f_clearcoat = Spectrum(clearcoat, clearcoat, clearcoat);

//             // Spectrum inital = Spectrum(1, 1, 1);
//             // Real F_c, D_c, G_c;
//             // Real R0 = Real(0.04);
//             // Real alpha_g = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);
//             // Real alpha_sq = alpha_g * alpha_g;
//             // F_c = R0 + (1 - R0) * pow(1 - abs(h_dot_out), Real(5));
//             // D_c = (alpha_sq - 1) / (c_PI * log(alpha_sq) * (1 + (alpha_sq - 1) * h_local.z * h_local.z));
//             // Real likeA_in, likeA_out;
//             // likeA_in = (sqrt(1 + (pow(local_dir_in.x * Real(0.25), Real(2)) + pow(local_dir_in.y * Real(0.25), Real(2))) / pow(local_dir_in.z, Real(2))) - 1) * Real(0.5);
//             // likeA_out = (sqrt(1 + (pow(local_dir_out.x * Real(0.25), Real(2)) + pow(local_dir_out.y * Real(0.25), Real(2))) / pow(local_dir_out.z, Real(2))) - 1) * Real(0.5);
//             // G_c = (1 / (1 + likeA_in)) * (1 / (1 + likeA_out));
//             // f_clearcoat = inital * F_c * D_c * G_c / (4 * abs(n_dot_in));
//         }
//     }

//     // glass calculation
//     {
//         Real F = fresnel_dielectric(h_dot_in, eta);
//         Real D = GTR2(h_local, roughness, anisotropic);
//         Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
//                  smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);
//         if (reflect) {
//             f_glass = base_color * (F * D * G) / (4 * abs(n_dot_in));
//         } else {
//             Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
//             Real sqrt_denom = h_dot_in + eta * h_dot_out;
//             f_glass = sqrt(base_color) * (eta_factor * (1 - F) * D * G * eta * eta * abs(h_dot_out * h_dot_in)) /
//                       (abs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom);
//         }
//     }
//     return diffuseWeight * f_diffuse + metalWeight * f_metal + glassWeight * f_glass + clearcoatWeight * f_clearcoat + sheen * (1 - metallic) * f_sheen;
// }

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                       dot(vertex.geometric_normal, dir_out) >
                   0;
    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;

    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    specular_transmission = std::clamp(specular_transmission, Real(0.01), Real(1));
    metallic = std::clamp(metallic, Real(0.01), Real(1));
    subsurface = std::clamp(subsurface, Real(0.01), Real(1));
    specular = std::clamp(specular, Real(0.01), Real(1));
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    specular_tint = std::clamp(specular_tint, Real(0.01), Real(1));
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    sheen = std::clamp(sheen, Real(0.01), Real(1));
    sheen_tint = std::clamp(sheen_tint, Real(0.01), Real(1));
    clearcoat = std::clamp(clearcoat, Real(0.01), Real(1));
    clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));

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
    Real n_dot_out = dot(frame.n, dir_out);

    Spectrum f_diffuse, f_sheen, f_metal, f_clearcoat, f_glass;

    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_sheen = (1 - metallic) * sheen;
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;

    auto R_0 = [](Real eta) {
        return ((eta - 1) * (eta - 1)) / ((eta + 1) * (eta + 1));
    };

    DisneyGlass bsdf_glass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
    f_glass = eval_op::operator()(bsdf_glass);

    if (inside) {
        return w_glass * f_glass;
    }

    DisneyDiffuse bsdf_diffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    f_diffuse = eval_op::operator()(bsdf_diffuse);

    DisneySheen bsdf_sheen{bsdf.base_color, bsdf.sheen_tint};
    f_sheen = eval_op::operator()(bsdf_sheen);

    DisneyClearcoat bsdf_clearcoat{bsdf.clearcoat_gloss};
    f_clearcoat = eval_op::operator()(bsdf_clearcoat);

    {
        // metal
        Spectrum C_tint = luminance(base_color) > 0 ? base_color / luminance(base_color) : Spectrum(1, 1, 1);
        Spectrum K_s = (1 - specular_tint) * Spectrum(1, 1, 1) + specular_tint * C_tint;
        Spectrum C_0 = specular * R_0(eta) * (1 - metallic) * K_s + metallic * base_color;
        // Spectrum F = schlick_fresnel(C_0, abs(h_dot_out));
        Spectrum F = C_0 + (1 - C_0) * pow(1 - abs(h_dot_out), 5);

        Real D = GTR2(h_local, roughness, anisotropic);

        Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
                 smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);

        f_metal = (F * D * G) / (4 * abs(n_dot_in));
    }

    return w_diffuse * f_diffuse +
           w_sheen * f_sheen +
           w_metal * f_metal +
           w_clearcoat * f_clearcoat +
           w_glass * f_glass;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                       dot(vertex.geometric_normal, dir_out) >
                   0;
    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    specular_transmission = std::clamp(specular_transmission, Real(0.01), Real(1));
    metallic = std::clamp(metallic, Real(0.01), Real(1));
    subsurface = std::clamp(subsurface, Real(0.01), Real(1));
    specular = std::clamp(specular, Real(0.01), Real(1));
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    specular_tint = std::clamp(specular_tint, Real(0.01), Real(1));
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    sheen = std::clamp(sheen, Real(0.01), Real(1));
    sheen_tint = std::clamp(sheen_tint, Real(0.01), Real(1));
    clearcoat = std::clamp(clearcoat, Real(0.01), Real(1));
    clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));

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
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, h);

    Real pdf_diffuse, pdf_metal, pdf_clearcoat, pdf_glass;

    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;
    Real w_sum = w_diffuse + w_metal + w_clearcoat + w_glass;
    w_diffuse /= w_sum;
    w_metal /= w_sum;
    w_clearcoat /= w_sum;
    w_glass /= w_sum;

    {
        // glass
        Real F = fresnel_dielectric(h_dot_in, eta);
        Real D = GTR2(h_local, roughness, anisotropic);
        Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic);

        if (reflect) {
            pdf_glass = (F * D * G) / (4 * abs(n_dot_in));
        } else {
            Real sqrt_denom = h_dot_in + eta * h_dot_out;
            Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
            pdf_glass = (1 - F) * D * G * abs(dh_dout * h_dot_in / n_dot_in);
        }
    }

    if (inside) {
        return pdf_glass;
        // return w_glass * pdf_glass;
    }

    {
        // diffuse
        // pdf_diffuse = fmax(dot(frame.n, dir_out), Real(0)) / c_PI;

        if (dot(vertex.geometric_normal, dir_out) > 0) {
            // Diffuse calculation
            {
                pdf_diffuse = fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
            }

        } else {
            pdf_diffuse = 0;
        }
    }
    {
        // metal
        Real D = GTR2(h_local, roughness, anisotropic);
        Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic);
        pdf_metal = (D * G) / (4 * abs(n_dot_in));
    }
    {
        // clearcoat
        Real alpha = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
        Real D = (alpha * alpha - 1) / (c_PI * log(alpha * alpha) * (1 + (alpha * alpha - 1) * h_local.z * h_local.z));

        pdf_clearcoat = (D * abs(n_dot_h)) / (4 * abs(h_dot_out));
    }

    return w_diffuse * pdf_diffuse +
           w_metal * pdf_metal +
           w_clearcoat * pdf_clearcoat +
           w_glass * pdf_glass;
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    specular_transmission = std::clamp(specular_transmission, Real(0.01), Real(1));
    metallic = std::clamp(metallic, Real(0.01), Real(1));
    subsurface = std::clamp(subsurface, Real(0.01), Real(1));
    specular = std::clamp(specular, Real(0.01), Real(1));
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    specular_tint = std::clamp(specular_tint, Real(0.01), Real(1));
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    sheen = std::clamp(sheen, Real(0.01), Real(1));
    sheen_tint = std::clamp(sheen_tint, Real(0.01), Real(1));
    clearcoat = std::clamp(clearcoat, Real(0.01), Real(1));
    clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 local_micro_normal = sample_visible_normals(to_local(frame, dir_in), alpha_x, alpha_y, rnd_param_uv);
    Vector3 h = to_world(frame, local_micro_normal);

    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;
    Real w_sum = w_diffuse + w_metal + w_clearcoat + w_glass;
    Real P_glass = w_glass / w_sum;
    Real P_diffuse = (w_glass + w_diffuse) / w_sum;
    Real P_metal = (w_glass + w_diffuse + w_metal) / w_sum;
    Real P_clearcoat = (w_glass + w_diffuse + w_metal + w_clearcoat) / w_sum;

    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    if (inside) {
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
            Vector3 refracted = -dir_in / eta + (abs(h_dot_in) / eta - h_dot_out) * h;
            return BSDFSampleRecord{refracted, eta, roughness};
        }
    } else if (rnd_param_w <= P_glass) {
        // Now we need to decide whether to reflect or refract.
        // We do this using the Fresnel term.
        Real h_dot_in = dot(h, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);

        if (rnd_param_w / P_glass <= F) {
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
            Vector3 refracted = -dir_in / eta + (abs(h_dot_in) / eta - h_dot_out) * h;
            return BSDFSampleRecord{refracted, eta, roughness};
        }
    } else if (rnd_param_w <= P_diffuse) {
        return BSDFSampleRecord{
            to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
            Real(0) /* eta */, eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool) /* roughness */};
    } else if (rnd_param_w <= P_metal) {
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);
        return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, roughness /* roughness */
        };
    } else if (rnd_param_w <= P_clearcoat) {
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
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
