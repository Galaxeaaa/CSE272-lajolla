#pragma once

#include "scene.h"
#include "pcg.h"

void updateMedium(const PathVertex &vertex, Vector3 dir, int &medium) {
    if (vertex.interior_medium_id != vertex.exterior_medium_id) {
        bool getting_out = dot(vertex.geometric_normal, dir) > 0;
        medium = getting_out ? vertex.exterior_medium_id : vertex.interior_medium_id;
    }
}

std::pair<const Light &, int> sampleLight(const Scene &scene, pcg32_state &rng) {
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    return {scene.lights[light_id], light_id};
}

PointAndNormal samplePointOnLight(const Scene &scene, const Light &light, const Vector3 &ref_point, pcg32_state &rng) {
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real shape_w = next_pcg32_real<Real>(rng);
    PointAndNormal p_light = sample_point_on_light(light, ref_point, light_uv, shape_w, scene);
    return p_light;
}

// The simplest volumetric renderer:
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);

    std::optional<PathVertex> vertex_ = intersect(scene, ray);
    if (!vertex_) {
        return make_zero_spectrum();
    }

    PathVertex vertex = *vertex_;
    Real t = distance(ray.org, vertex.position);
    Spectrum sigma_a = get_sigma_a(scene.media[0], vertex.position);
    Spectrum transmittance = exp(-sigma_a * t);
    Spectrum Le = make_zero_spectrum();

    if (is_light(scene.shapes[vertex.shape_id])) {
        Le = emission(vertex, -ray.dir, scene);
    }

    return transmittance * Le;
}

// The second simplest volumetric renderer:
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    Real t_hit = infinity<Real>();
    PathVertex vertex;
    if (vertex_) {
        vertex = *vertex_;
        t_hit = distance(ray.org, vertex.position);
    }

    Spectrum sigma_a = get_sigma_a(scene.media[0], Vector3(0, 0, 0));
    Spectrum sigma_s = get_sigma_s(scene.media[0], Vector3(0, 0, 0));
    Spectrum sigma_t_v = sigma_a + sigma_s;
    Real sigma_t = sigma_a.x + sigma_s.x;  // the volume is monochromatic

    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1 - u) / sigma_t;
    Spectrum transmittance;

    if (t < t_hit) {
        // if (true) {
        // Monte Carlo sampling transmittance
        Real pdf_transmittance = exp(-sigma_t * t) * sigma_t;
        transmittance = exp(-sigma_t_v * t);

        // current point p
        Vector3 p = ray.org + t * ray.dir;

        auto [light, light_id] = sampleLight(scene, rng);
        PointAndNormal p_light = samplePointOnLight(scene, light, p, rng);

        Vector3 dir_out = normalize(p_light.position - p);
        Spectrum phase = eval(get_phase_function(scene.media[0]), -ray.dir, dir_out);  // only 1 volume
        Spectrum Le_light = emission(light,
                                     -dir_out,
                                     Real(0),
                                     p_light,
                                     scene);

        Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p, scene);

        Spectrum L_s1_estimate = phase * Le_light * exp(-sigma_t * distance(p, p_light.position)) * (abs(dot(dir_out, p_light.normal)) / distance_squared(p, p_light.position));

        return (transmittance / pdf_transmittance) * sigma_s * (L_s1_estimate / L_s1_pdf);
    } else {
        // Real pdf_transmittance = exp(-sigma_t * t_hit) + sigma_t;
        Real pdf_transmittance = exp(-sigma_t * t_hit);
        transmittance = exp(-sigma_t_v * t_hit);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        return (transmittance / pdf_transmittance) * Le;
    }

    return make_zero_spectrum();
}

Spectrum vol_path_tracing_2_bonus(const Scene &scene,
                                  int x, int y, /* pixel coordinates */
                                  pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    Real t_hit = infinity<Real>();
    PathVertex vertex;
    if (vertex_) {
        vertex = *vertex_;
        t_hit = distance(ray.org, vertex.position);
    }

    Spectrum sigma_a = get_sigma_a(scene.media[0], Vector3(0, 0, 0));
    Spectrum sigma_s = get_sigma_s(scene.media[0], Vector3(0, 0, 0));
    Spectrum sigma_t_v = sigma_a + sigma_s;
    Real sigma_t = sigma_a.x + sigma_s.x;  // the volume is monochromatic

    Spectrum radiance = make_zero_spectrum();

    {
        auto [light, light_id] = sampleLight(scene, rng);

        Real D = 0;
        int shape_id = std::get<DiffuseAreaLight>(light).shape_id;
        Vector3 p_light_center = std::get<Sphere>(scene.shapes[shape_id]).position;
        Vector3 v_light = p_light_center - ray.org;
        Real Delta = dot(v_light, ray.dir);
        D = sqrt(length_squared(v_light) - Delta * Delta);
        Real theta_a = atan(-Delta / D);
        Real theta_b = atan((infinity<Real>() - Delta) / D);
        Real u = next_pcg32_real<Real>(rng);
        Real t = D * tan((1 - u) * theta_a + u * theta_b);
        Vector3 p = ray.org + (Delta + t) * ray.dir;

        if (Delta + t < t_hit) {
            auto [light, light_id] = sampleLight(scene, rng);
            PointAndNormal p_light = samplePointOnLight(scene, light, p, rng);

            Vector3 dir_out = normalize(p_light.position - p);
            Spectrum phase = eval(get_phase_function(scene.media[0]), -ray.dir, dir_out);  // only 1 volume
            Spectrum Le_light = emission(light, -dir_out, Real(0), p_light, scene);

            Spectrum transmittance = exp(-sigma_t_v * (Delta + t));
            Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p, scene);
            Spectrum L_s1_estimate = phase * Le_light * exp(-sigma_t * distance(p, p_light.position)) * (abs(dot(dir_out, p_light.normal)) / distance_squared(p, p_light.position));

            Real pdf_equi = D / (theta_b - theta_a) / (D * D + t * t);
            Real pdf_trans = exp(-sigma_t * (Delta + t)) * sigma_t;
            Real w = pdf_equi * pdf_equi / (pdf_equi * pdf_equi + pdf_trans * pdf_trans);

            radiance += transmittance * sigma_s * (L_s1_estimate / L_s1_pdf) / pdf_equi * w;
        } else {
            Spectrum transmittance = exp(-sigma_t_v * t_hit);
            Real pdf_equi = (atan(infinity<Real>() / D) - atan((t_hit - Delta) / D)) / (theta_b - theta_a);
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }

            Real pdf_trans = exp(-sigma_t * t_hit);
            Real w = pdf_equi * pdf_equi / (pdf_equi * pdf_equi + pdf_trans * pdf_trans);

            radiance += (transmittance / pdf_equi) * Le * w;
        }
    }

    {
        Real u = next_pcg32_real<Real>(rng);
        Real t = -log(1 - u) / sigma_t;
        Spectrum transmittance;

        auto [light, light_id] = sampleLight(scene, rng);
        int shape_id = std::get<DiffuseAreaLight>(light).shape_id;
        Vector3 p_light_center = std::get<Sphere>(scene.shapes[shape_id]).position;

        if (t < t_hit) {
            Real pdf_transmittance = exp(-sigma_t * t) * sigma_t;
            transmittance = exp(-sigma_t_v * t);
            Vector3 p = ray.org + t * ray.dir;

            auto [light, light_id] = sampleLight(scene, rng);
            PointAndNormal p_light = samplePointOnLight(scene, light, p, rng);

            Vector3 dir_out = normalize(p_light.position - p);
            Spectrum phase = eval(get_phase_function(scene.media[0]), -ray.dir, dir_out);  // only 1 volume
            Spectrum Le_light = emission(light,
                                         -dir_out,
                                         Real(0),
                                         p_light,
                                         scene);

            Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p, scene);

            Spectrum L_s1_estimate = phase * Le_light * exp(-sigma_t * distance(p, p_light.position)) * (abs(dot(dir_out, p_light.normal)) / distance_squared(p, p_light.position));

            Real D = 0;
            Vector3 v_light = p_light_center - ray.org;
            Real Delta = dot(v_light, ray.dir);
            D = sqrt(length_squared(v_light) - Delta * Delta);
            Real theta_a = atan(-Delta / D);
            Real theta_b = atan((infinity<Real>() - Delta) / D);
            Real pdf_equi = D / (theta_b - theta_a) / (D * D + (t - Delta) * (t - Delta));
            Real w = pdf_transmittance * pdf_transmittance / (pdf_transmittance * pdf_transmittance + pdf_equi * pdf_equi);

            radiance += (transmittance / pdf_transmittance) * sigma_s * (L_s1_estimate / L_s1_pdf) * w;
        } else {
            // Real pdf_transmittance = exp(-sigma_t * t_hit) + sigma_t;
            Real pdf_transmittance = exp(-sigma_t * t_hit);
            transmittance = exp(-sigma_t_v * t_hit);
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }

            Real D = 0;
            Vector3 v_light = p_light_center - ray.org;
            Real Delta = dot(v_light, ray.dir);
            D = sqrt(length_squared(v_light) - Delta * Delta);
            Real theta_a = atan(-Delta / D);
            Real theta_b = atan((infinity<Real>() - Delta) / D);
            Real pdf_equi = (atan((infinity<Real>() - Delta) / D) - atan((t_hit - Delta) / D)) / (theta_b - theta_a);
            Real w = pdf_transmittance * pdf_transmittance / (pdf_transmittance * pdf_transmittance + pdf_equi * pdf_equi);

            radiance += (transmittance / pdf_transmittance) * Le * w;
        }
    }

    return radiance;
}

// The third volumetric renderer (not so simple anymore):
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium = scene.camera.medium_id;

    Spectrum current_path_throughput = Spectrum(1, 1, 1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = infinity<Real>();
        PathVertex vertex;
        if (vertex_) {
            vertex = *vertex_;
            t_hit = distance(ray.org, vertex.position);
        }
        Spectrum transmittance = Spectrum(1, 1, 1);
        Real pdf_transmittance = 1;

        if (current_medium >= 0) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t_v = sigma_a + sigma_s;
            Real sigma_t = sigma_a.x + sigma_s.x;  // the volume is monochromatic
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            if (t < t_hit) {
                pdf_transmittance = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t_v * t);
                scatter = true;
                ray.org = ray.org + t * ray.dir;
            } else {
                pdf_transmittance = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t_v * t_hit);
            }
        }

        current_path_throughput *= transmittance / pdf_transmittance;

        // hit object
        if (!scatter && vertex.shape_id >= 0) {
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }
            radiance += current_path_throughput * Le;
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            break;
        }

        // hit index-matching surface
        if (!scatter && vertex_) {
            if (vertex.material_id == -1) {
                if (vertex.interior_medium_id != vertex.exterior_medium_id) {
                    bool getting_out = dot(vertex.geometric_normal, ray.dir) > 0;
                    current_medium = getting_out ? vertex.exterior_medium_id : vertex.interior_medium_id;
                }
                ray.org = vertex.position + ray.dir * get_intersection_epsilon(scene);
                bounces++;
                continue;
            }
        }

        if (scatter) {
            Vector2 phase_uv = Vector2{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);
            auto next_dir_ = sample_phase_function(phase_f, -ray.dir, phase_uv);
            Vector3 next_dir;
            if (next_dir_) {
                next_dir = *next_dir_;
            }
            Spectrum phase = eval(phase_f, -ray.dir, next_dir);
            Real pdf_phase = pdf_sample_phase(phase_f, -ray.dir, next_dir);
            current_path_throughput *= (phase / pdf_phase) * get_sigma_s(scene.media[current_medium], ray.org);
            ray.dir = next_dir;
        } else {
            // hit a surface
            break;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(length(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    return radiance;
}

// The fourth volumetric renderer:
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    int current_medium = scene.camera.medium_id;

    Spectrum current_path_throughput = Spectrum(1, 1, 1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    // Caches for calculating pdf of NEE.
    Real pdf_sampling_phase_cache = 0;
    Vector3 pos_cache = ray.org;  // for NEE
    Spectrum phase_cache = make_zero_spectrum();
    Real pdf_multi_trans_cache = 1;
    bool never_scatter = true;

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray);
        Real t_hit = infinity<Real>();
        PathVertex vertex;
        if (vertex_) {
            vertex = *vertex_;
            t_hit = distance(ray.org, vertex.position);
        }
        Spectrum transmittance = Spectrum(1, 1, 1);
        Real pdf_transmittance = 1;

        if (current_medium >= 0) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t_v = sigma_a + sigma_s;
            Real sigma_t = sigma_a.x + sigma_s.x;  // the volume is monochromatic
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            if (t < t_hit) {
                pdf_transmittance = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t_v * t);
                scatter = true;
                never_scatter = false;
                ray.org = ray.org + t * ray.dir;
            } else {
                pdf_transmittance = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t_v * t_hit);
                ray.org = ray.org + t_hit * ray.dir;
            }
            pdf_multi_trans_cache *= pdf_transmittance;
        }

        current_path_throughput *= transmittance / pdf_transmittance;

        // hit object
        if (!scatter && vertex.material_id >= 0) {
            if (never_scatter) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    Spectrum Le = emission(vertex, -ray.dir, scene);
                    radiance += current_path_throughput * Le;
                }
            } else {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    int id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Real pdf_nee = light_pmf(scene, id) * pdf_point_on_light(scene.lights[id],
                                                                             PointAndNormal{vertex.position, vertex.geometric_normal},
                                                                             pos_cache,
                                                                             scene);
                    Vector3 dir_out = normalize(vertex.position - pos_cache);
                    Real G = abs(dot(dir_out, vertex.geometric_normal)) / distance_squared(pos_cache, vertex.position);
                    Real pdf_dir = pdf_sampling_phase_cache * pdf_multi_trans_cache * G;
                    Real w = (pdf_dir * pdf_dir) / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                    // Real w = 0.0;
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            break;
        }

        // hit index-matching surface
        if (!scatter && vertex_) {
            if (vertex.material_id == -1) {
                if (vertex.interior_medium_id != vertex.exterior_medium_id) {
                    bool getting_out = dot(vertex.geometric_normal, ray.dir) > 0;
                    current_medium = getting_out ? vertex.exterior_medium_id : vertex.interior_medium_id;
                }
                ray.org = vertex.position + ray.dir * get_intersection_epsilon(scene);
                bounces++;
                continue;
            }
        }

        if (scatter) {
            // NEE
            auto [light, light_id] = sampleLight(scene, rng);
            PointAndNormal p_light = samplePointOnLight(scene, light, ray.org, rng);

            Vector3 p_org = ray.org;
            Vector3 p = p_org;
            Real T_light = 1;
            int shadow_medium = current_medium;
            int shadow_bounces = 0;
            Real p_trans_dir = 1;
            while (true) {
                Ray shadow_ray = Ray{p, normalize(p_light.position - p), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_light.position, p)};
                // Ray shadow_ray = Ray{p, normalize(p_light.position - p), Real(0), infinity<Real>()};
                std::optional<PathVertex> intersection_ = intersect(scene, shadow_ray);
                Real t = distance(p, p_light.position);
                if (intersection_) {
                    t = distance(p, intersection_->position);
                }
                if (shadow_medium >= 0) {
                    Real sigma_t = get_sigma_a(scene.media[shadow_medium], p).x + get_sigma_s(scene.media[shadow_medium], p).x;
                    T_light *= exp(-sigma_t * t);
                    p_trans_dir *= exp(-sigma_t * t);
                }
                if (!intersection_) {
                    break;
                } else {
                    if (intersection_->material_id >= 0) {
                        // return make_zero_spectrum();
                        T_light = 0;
                        break;
                    }
                    shadow_bounces++;
                    if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                        // return make_zero_spectrum();
                        T_light = 0;
                        break;
                    }
                    if (intersection_->material_id == -1) {
                        if (intersection_->interior_medium_id != intersection_->exterior_medium_id) {
                            bool getting_out = dot(intersection_->geometric_normal, shadow_ray.dir) > 0;
                            shadow_medium = getting_out ? intersection_->exterior_medium_id : intersection_->interior_medium_id;
                        }
                    }
                    p = p + t * shadow_ray.dir;
                }
            }

            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], p_org);
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);

            if (T_light > 0) {
                Vector3 dir_out = normalize(p_light.position - p_org);
                Real G = abs(dot(dir_out, p_light.normal)) / distance_squared(p_org, p_light.position);
                Spectrum Le_light = emission(light, -dir_out, Real(0), p_light, scene);
                Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p_org, scene);
                Spectrum phase = eval(phase_f, -ray.dir, dir_out);
                Spectrum contribution = T_light * G * phase * Le_light / pdf_nee;
                Real pdf_phase = pdf_sample_phase(phase_f, -ray.dir, dir_out) * G * p_trans_dir;
                Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
                radiance += current_path_throughput * sigma_s * contribution * w;
            }

            // phase sampling
            Vector2 phase_uv = Vector2{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto next_dir_ = sample_phase_function(phase_f, -ray.dir, phase_uv);
            Vector3 next_dir;
            if (next_dir_) {
                next_dir = *next_dir_;
            }
            Spectrum phase = eval(phase_f, -ray.dir, next_dir);
            Real pdf_phase = pdf_sample_phase(phase_f, -ray.dir, next_dir);

            pos_cache = ray.org;
            phase_cache = phase;
            pdf_sampling_phase_cache = pdf_phase;
            pdf_multi_trans_cache = 1;
            current_path_throughput *= (phase / pdf_phase) * sigma_s;
            ray.dir = next_dir;
        } else {
            // hit a surface
            break;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(length(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    return radiance;
}

// The fifth volumetric renderer:
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    int current_medium = scene.camera.medium_id;

    Spectrum current_path_throughput = Spectrum(1, 1, 1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    // Caches for calculating pdf of NEE.
    Vector3 pos_cache = ray.org;  // for NEE
    Real pdf_sampling_cache = 0;
    Real pdf_multi_trans_cache = 1;
    bool never_scatter_or_bsdf = true;

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray);
        Real t_hit = infinity<Real>();
        PathVertex vertex;
        if (vertex_) {
            vertex = *vertex_;
            t_hit = distance(ray.org, vertex.position);
        }
        Spectrum transmittance = make_const_spectrum(Real(1));
        Real pdf_transmittance = 1;

        if (current_medium >= 0) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t_v = sigma_a + sigma_s;
            Real sigma_t = sigma_a.x + sigma_s.x;  // the volume is monochromatic
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            if (t < t_hit) {
                pdf_transmittance = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t_v * t);
                scatter = true;
                ray.org = ray.org + t * ray.dir;
            } else {
                pdf_transmittance = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t_v * t_hit);
                ray.org = ray.org + t_hit * ray.dir;
            }
            pdf_multi_trans_cache *= pdf_transmittance;
        }

        current_path_throughput *= transmittance / pdf_transmittance;

        // hit object
        if (!scatter && vertex_) {
            if (is_light(scene.shapes[vertex.shape_id])) {
                if (never_scatter_or_bsdf) {
                    Spectrum Le = emission(vertex, -ray.dir, scene);
                    radiance += current_path_throughput * Le;
                } else {
                    // need to consider MIS
                    int id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Real pdf_nee = light_pmf(scene, id) * pdf_point_on_light(scene.lights[id],
                                                                             PointAndNormal{vertex.position, vertex.geometric_normal},
                                                                             pos_cache,
                                                                             scene);
                    Vector3 dir_out = normalize(vertex.position - pos_cache);
                    Real G = abs(dot(dir_out, vertex.geometric_normal)) / distance_squared(pos_cache, vertex.position);
                    Real pdf_dir = pdf_sampling_cache * pdf_multi_trans_cache * G;
                    Real w = (pdf_dir * pdf_dir) / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            break;
        }

        if (!scatter && vertex_) {
            if (vertex.material_id == -1) {
                // Hit index-matching surface
                // ---------------------------------
                updateMedium(vertex, ray.dir, current_medium);
                Vector3 N = dot(ray.dir, vertex.geometric_normal) > 0 ? vertex.geometric_normal : -vertex.geometric_normal;
                ray.org = vertex.position + N * get_intersection_epsilon(scene);
                bounces++;
                continue;
            } else {
                // Hit surface
                // ---------------------------------

                // NEE
                Vector3 p_org = vertex.position;
                auto [light, light_id] = sampleLight(scene, rng);
                PointAndNormal p_light = samplePointOnLight(scene, light, p_org, rng);

                Vector3 p = p_org;
                Vector3 dir_in = -ray.dir;
                Vector3 dir_out = normalize(p_light.position - p_org);
                Real T_light = 1;
                Real p_trans_dir = 1;
                int shadow_medium = current_medium;
                int shadow_bounces = 0;
                while (true) {
                    Ray shadow_ray = Ray{p, dir_out, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_light.position, p)};
                    std::optional<PathVertex> intersection_ = intersect(scene, shadow_ray);
                    Real t = distance(p, p_light.position);
                    if (intersection_) {
                        t = distance(p, intersection_->position);
                    }
                    if (shadow_medium >= 0) {
                        Real sigma_t = get_sigma_a(scene.media[shadow_medium], p).x + get_sigma_s(scene.media[shadow_medium], p).x;
                        T_light *= exp(-sigma_t * t);
                        p_trans_dir *= exp(-sigma_t * t);
                    }
                    if (!intersection_) {
                        break;
                    } else {
                        PathVertex intersection = *intersection_;
                        if (intersection.material_id >= 0) {
                            // hit surface
                            T_light = 0;
                            break;
                        }
                        shadow_bounces++;
                        if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                            T_light = 0;
                            break;
                        }
                        if (intersection.material_id == -1) {
                            // hit index-matching surface
                            updateMedium(intersection, shadow_ray.dir, shadow_medium);
                        }
                        p = p + t * shadow_ray.dir;
                    }
                }

                const Material &mat = scene.materials[vertex.material_id];

                if (T_light > 0) {
                    Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p_org, scene);
                    if (pdf_nee <= 0) {
                        break;
                    }
                    Vector3 dir_out = normalize(p_light.position - p_org);
                    Real G = abs(dot(dir_out, p_light.normal)) / distance_squared(p_org, p_light.position);
                    Spectrum Le_light = emission(light, -dir_out, Real(0), p_light, scene);
                    Spectrum f = eval(mat, dir_in, dir_out, vertex, scene.texture_pool);
                    Real pdf_bsdf = G * pdf_sample_bsdf(mat, dir_in, dir_out, vertex, scene.texture_pool);
                    Spectrum contribution = T_light * G * f * Le_light / pdf_nee;
                    Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
                    radiance += current_path_throughput * contribution * w;
                }

                // BSDF sampling
                Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_ = sample_bsdf(mat, dir_in, vertex, scene.texture_pool, bsdf_rnd_param_uv, bsdf_rnd_param_w);
                if (!bsdf_sample_) {
                    break;
                }
                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
                Vector3 next_dir = bsdf_sample.dir_out;
                Spectrum bsdf = eval(mat, dir_in, next_dir, vertex, scene.texture_pool);
                Real pdf_bsdf = pdf_sample_bsdf(mat, dir_in, next_dir, vertex, scene.texture_pool);
                if (pdf_bsdf <= 0) {
                    break;
                }
                current_path_throughput *= bsdf / pdf_bsdf;
                pdf_sampling_cache = pdf_bsdf;

                pdf_multi_trans_cache = 1;
                Vector3 N = dot(next_dir, vertex.geometric_normal) > 0 ? vertex.geometric_normal : -vertex.geometric_normal;
                ray.org = vertex.position + N * get_intersection_epsilon(scene);
                ray.dir = next_dir;
                updateMedium(vertex, ray.dir, current_medium);
                pos_cache = ray.org;
                never_scatter_or_bsdf = false;
            }
        }

        if (scatter) {
            // Volume scatter
            // --------------------------------

            // NEE
            auto [light, light_id] = sampleLight(scene, rng);
            PointAndNormal p_light = samplePointOnLight(scene, light, ray.org, rng);

            Vector3 p = ray.org;
            Real T_light = 1;
            int shadow_medium = current_medium;
            int shadow_bounces = 0;
            Real p_trans_dir = 1;
            while (true) {
                Ray shadow_ray = Ray{p, normalize(p_light.position - p), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_light.position, p)};
                std::optional<PathVertex> intersection_ = intersect(scene, shadow_ray);
                Real t = distance(p, p_light.position);
                if (intersection_) {
                    t = distance(p, intersection_->position);
                }
                if (shadow_medium >= 0) {
                    Real sigma_t = get_sigma_a(scene.media[shadow_medium], p).x + get_sigma_s(scene.media[shadow_medium], p).x;
                    T_light *= exp(-sigma_t * t);
                    p_trans_dir *= exp(-sigma_t * t);
                }
                if (!intersection_) {
                    break;
                } else {
                    PathVertex intersection = *intersection_;
                    if (intersection.material_id >= 0) {
                        T_light = 0;
                        break;
                    }
                    shadow_bounces++;
                    if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                        T_light = 0;
                        break;
                    }
                    if (intersection.material_id == -1) {
                        updateMedium(intersection, shadow_ray.dir, shadow_medium);
                    }
                    p = p + t * shadow_ray.dir;
                }
            }

            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);

            if (T_light > 0) {
                Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, ray.org, scene);
                if (pdf_nee <= 0) {
                    break;
                }
                Vector3 dir_out = normalize(p_light.position - ray.org);
                Real G = abs(dot(dir_out, p_light.normal)) / distance_squared(ray.org, p_light.position);
                Spectrum Le_light = emission(light, -dir_out, Real(0), p_light, scene);
                Spectrum phase = eval(phase_f, -ray.dir, dir_out);
                Spectrum contribution = T_light * G * phase * Le_light / pdf_nee;
                Real pdf_phase = pdf_sample_phase(phase_f, -ray.dir, dir_out) * G * p_trans_dir;
                Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
                radiance += current_path_throughput * sigma_s * contribution * w;
            }

            // Phase sampling
            Vector2 phase_uv = Vector2{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto next_dir_ = sample_phase_function(phase_f, -ray.dir, phase_uv);
            Vector3 next_dir;
            if (next_dir_) {
                next_dir = *next_dir_;
            }
            Spectrum phase = eval(phase_f, -ray.dir, next_dir);
            Real pdf_phase = pdf_sample_phase(phase_f, -ray.dir, next_dir);

            current_path_throughput *= (phase / pdf_phase) * sigma_s;
            pos_cache = ray.org;
            pdf_sampling_cache = pdf_phase;
            pdf_multi_trans_cache = 1;
            ray.dir = next_dir;
            never_scatter_or_bsdf = false;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(length(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    return radiance;
}

// The final volumetric renderer:
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    int current_medium = scene.camera.medium_id;

    Spectrum current_path_throughput = make_const_spectrum(Real(1));
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    // Caches for calculating pdf of NEE.
    Vector3 pos_cache = ray.org;
    Real pdf_sampling_cache = 0;
    Spectrum pdf_trans_dir_cache = make_const_spectrum(Real(1));
    Spectrum pdf_trans_nee_cache = make_const_spectrum(Real(1));
    bool never_scatter_or_bsdf = true;

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray);
        Real t_hit = infinity<Real>();
        PathVertex vertex;
        if (vertex_) {
            vertex = *vertex_;
            t_hit = distance(ray.org, vertex.position);
        }
        Spectrum transmittance = make_const_spectrum(Real(1));
        Spectrum pdf_transmittance_dir = make_const_spectrum(Real(1));
        Spectrum pdf_transmittance_nee = make_const_spectrum(Real(1));

        if (current_medium >= 0) {
            const Medium &med = scene.media[current_medium];
            Spectrum sigma_m = get_majorant(med, ray);
            Real u = next_pcg32_real<Real>(rng);
            Real channel = std::clamp(int(u * 3), 0, 2);
            Real accum_t = 0;
            int iteration = 0;
            while (true) {
                if (sigma_m[channel] <= 0) {
                    break;
                }
                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }
                Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_m[channel];
                Real dt = t_hit - accum_t;
                accum_t = min(accum_t + t, t_hit);
                if (t < dt) {
                    Vector3 p = ray.org + accum_t * ray.dir;  // ?
                    Spectrum sigma_t = get_sigma_a(med, p) + get_sigma_s(med, p);
                    Spectrum real_prob = sigma_t / sigma_m;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        scatter = true;
                        transmittance *= exp(-sigma_m * t) / max(sigma_m);
                        pdf_transmittance_dir *= exp(-sigma_m * t) * sigma_m * real_prob / max(sigma_m);
                        // don't need to account for pdf_trans_nee since we scatter
                        ray.org = ray.org + accum_t * ray.dir;
                        break;
                    } else {
                        // hit a fake particle
                        transmittance *= exp(-sigma_m * t) * (sigma_m - sigma_t) / max(sigma_m);
                        pdf_transmittance_dir *= exp(-sigma_m * t) * sigma_m * (1 - real_prob) / max(sigma_m);
                        pdf_transmittance_nee *= exp(-sigma_m * t) * sigma_m / max(sigma_m);
                    }
                } else {
                    transmittance *= exp(-sigma_m * dt);
                    pdf_transmittance_dir *= exp(-sigma_m * dt);
                    pdf_transmittance_nee *= exp(-sigma_m * dt);
                    ray.org = ray.org + t_hit * ray.dir;
                    break;
                }
                iteration++;
            }

            pdf_trans_dir_cache *= pdf_transmittance_dir;
            pdf_trans_nee_cache *= pdf_transmittance_nee;
        }

        current_path_throughput *= transmittance / avg(pdf_transmittance_dir);

        // hit object
        if (!scatter && vertex_) {
            if (is_light(scene.shapes[vertex.shape_id])) {
                if (never_scatter_or_bsdf) {
                    Spectrum Le = emission(vertex, -ray.dir, scene);
                    radiance += current_path_throughput * Le;
                } else {
                    // need to consider MIS
                    int id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Real pdf_nee = light_pmf(scene, id) * pdf_point_on_light(scene.lights[id],
                                                                             PointAndNormal{vertex.position, vertex.geometric_normal},
                                                                             pos_cache,
                                                                             scene);
                    pdf_nee *= avg(pdf_trans_nee_cache);
                    Vector3 dir_out = normalize(vertex.position - pos_cache);
                    Real G = abs(dot(dir_out, vertex.geometric_normal)) / distance_squared(pos_cache, vertex.position);
                    Real pdf_dir = pdf_sampling_cache * avg(pdf_trans_dir_cache) * G;
                    Real w = (pdf_dir * pdf_dir) / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            break;
        }

        if (!scatter && vertex_) {
            if (vertex.material_id == -1) {
                // Hit index-matching surface
                // ---------------------------------
                updateMedium(vertex, ray.dir, current_medium);
                Vector3 N = dot(ray.dir, vertex.geometric_normal) > 0 ? vertex.geometric_normal : -vertex.geometric_normal;
                ray.org = vertex.position + N * get_intersection_epsilon(scene);
                bounces++;
                continue;
            } else {
                // Hit surface
                // ---------------------------------

                // NEE
                Vector3 p_org = vertex.position;
                auto [light, light_id] = sampleLight(scene, rng);
                PointAndNormal p_light = samplePointOnLight(scene, light, p_org, rng);

                Vector3 p = p_org;
                Vector3 dir_in = -ray.dir;
                Vector3 dir_out = normalize(p_light.position - p_org);
                Spectrum T_light = make_const_spectrum(1);
                Spectrum p_trans_dir = make_const_spectrum(1);
                Spectrum p_trans_nee = make_const_spectrum(1);
                int shadow_medium = current_medium;
                int shadow_bounces = 0;
                while (true) {
                    Ray shadow_ray = Ray{p, dir_out, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_light.position, p)};
                    std::optional<PathVertex> intersection_ = intersect(scene, shadow_ray);
                    Real next_t = distance(p, p_light.position);
                    if (intersection_) {
                        next_t = distance(p, intersection_->position);
                    }
                    if (shadow_medium >= 0) {
                        const Medium &med = scene.media[shadow_medium];
                        Spectrum sigma_m = get_majorant(med, shadow_ray);
                        Real u = next_pcg32_real<Real>(rng);
                        Real channel = std::clamp(int(u * 3), 0, 2);
                        Real accum_t = 0;
                        int iteration = 0;
                        while (true) {
                            if (sigma_m[channel] <= 0) {
                                break;
                            }
                            if (iteration >= scene.options.max_null_collisions) {
                                break;
                            }
                            Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_m[channel];
                            Real dt = next_t - accum_t;
                            accum_t = min(accum_t + t, next_t);
                            if (t < dt) {
                                // didn't hit surface, null-scattering
                                p = p + t * dir_out;
                                Spectrum sigma_t = get_sigma_a(med, p) + get_sigma_s(med, p);
                                T_light *= exp(-sigma_m * t) * (sigma_m - sigma_t) / max(sigma_m);
                                p_trans_nee *= exp(-sigma_m * t) * sigma_m / max(sigma_m);
                                Spectrum real_prob = sigma_t / sigma_m;
                                p_trans_dir *= exp(-sigma_m * t) * sigma_m * (1 - real_prob) / max(sigma_m);
                                if (max(T_light) <= 0) {
                                    break;
                                }
                            } else {
                                // reach next_t
                                p = p + dt * dir_out;
                                T_light *= exp(-sigma_m * dt);
                                p_trans_dir *= exp(-sigma_m * dt);
                                p_trans_nee *= exp(-sigma_m * dt);
                                break;
                            }
                            iteration++;
                        }
                    }
                    if (!intersection_) {
                        break;
                    } else {
                        PathVertex intersection = *intersection_;
                        if (intersection.material_id >= 0) {
                            T_light = make_zero_spectrum();
                            break;
                        }
                        shadow_bounces++;
                        if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                            T_light = make_zero_spectrum();
                            break;
                        }
                        if (intersection.material_id == -1) {
                            updateMedium(intersection, shadow_ray.dir, shadow_medium);
                        }
                        p = intersection.position;
                    }
                }

                const Material &mat = scene.materials[vertex.material_id];

                if (max(T_light) > 0) {
                    T_light.x = T_light.x > 0 ? T_light.x : 0;
                    T_light.y = T_light.y > 0 ? T_light.y : 0;
                    T_light.z = T_light.z > 0 ? T_light.z : 0;
                    Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p_org, scene);
                    if (pdf_nee <= 0) {
                        break;
                    }
                    pdf_nee *= avg(p_trans_nee);
                    Real G = abs(dot(dir_out, p_light.normal)) / distance_squared(p_org, p_light.position);
                    Spectrum Le_light = emission(light, -dir_out, Real(0), p_light, scene);
                    Spectrum f = eval(mat, dir_in, dir_out, vertex, scene.texture_pool);
                    Real pdf_bsdf = G * pdf_sample_bsdf(mat, dir_in, dir_out, vertex, scene.texture_pool) * avg(p_trans_dir);
                    Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
                    Spectrum contribution = T_light * G * f * Le_light / pdf_nee;
                    radiance += current_path_throughput * contribution * w;
                }

                // BSDF sampling
                Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_ = sample_bsdf(mat, dir_in, vertex, scene.texture_pool, bsdf_rnd_param_uv, bsdf_rnd_param_w);
                if (!bsdf_sample_) {
                    break;
                }
                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
                Vector3 next_dir = bsdf_sample.dir_out;
                Spectrum bsdf = eval(mat, dir_in, next_dir, vertex, scene.texture_pool);
                Real pdf_bsdf = pdf_sample_bsdf(mat, dir_in, next_dir, vertex, scene.texture_pool);
                if (pdf_bsdf <= 0) {
                    break;
                }
                current_path_throughput *= bsdf / pdf_bsdf;

                Vector3 N = dot(next_dir, vertex.geometric_normal) > 0 ? vertex.geometric_normal : -vertex.geometric_normal;
                ray.org = vertex.position + N * get_intersection_epsilon(scene);
                ray.dir = next_dir;
                updateMedium(vertex, ray.dir, current_medium);

                pos_cache = ray.org;
                pdf_sampling_cache = pdf_bsdf;
                pdf_trans_dir_cache = make_const_spectrum(1);
                pdf_trans_nee_cache = make_const_spectrum(1);
                never_scatter_or_bsdf = false;
            }
        }

        if (scatter) {
            // Volume scatter
            // --------------------------------
            // NEE
            Vector3 p_org = ray.org;
            auto [light, light_id] = sampleLight(scene, rng);
            PointAndNormal p_light = samplePointOnLight(scene, light, p_org, rng);

            Vector3 p = p_org;
            Vector3 dir_in = -ray.dir;
            Vector3 dir_out = normalize(p_light.position - p_org);
            Spectrum T_light = make_const_spectrum(1);
            Spectrum p_trans_dir = make_const_spectrum(1);
            Spectrum p_trans_nee = make_const_spectrum(1);
            int shadow_medium = current_medium;
            int shadow_bounces = 0;
            while (true) {
                Ray shadow_ray = Ray{p, dir_out, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_light.position, p)};
                std::optional<PathVertex> intersection_ = intersect(scene, shadow_ray);
                Real next_t = distance(p, p_light.position);
                if (intersection_) {
                    next_t = distance(p, intersection_->position);
                }
                if (shadow_medium >= 0) {
                    const Medium &med = scene.media[shadow_medium];
                    Spectrum sigma_m = get_majorant(med, shadow_ray);
                    Real u = next_pcg32_real<Real>(rng);
                    Real channel = std::clamp(int(u * 3), 0, 2);
                    Real accum_t = 0;
                    int iteration = 0;
                    while (true) {
                        if (sigma_m[channel] <= 0) {
                            break;
                        }
                        if (iteration >= scene.options.max_null_collisions) {
                            break;
                        }
                        Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_m[channel];
                        Real dt = next_t - accum_t;
                        accum_t = min(accum_t + t, next_t);
                        if (t < dt) {
                            // didn't hit surface, null-scattering
                            p = p + t * dir_out;
                            Spectrum sigma_t = get_sigma_a(med, p) + get_sigma_s(med, p);
                            T_light *= exp(-sigma_m * t) * (sigma_m - sigma_t) / max(sigma_m);
                            p_trans_nee *= exp(-sigma_m * t) * sigma_m / max(sigma_m);
                            Spectrum real_prob = sigma_t / sigma_m;
                            p_trans_dir *= exp(-sigma_m * t) * sigma_m * (1 - real_prob) / max(sigma_m);
                            if (max(T_light) <= 0) {
                                break;
                            }
                        } else {
                            // reach next_t
                            p = p + dt * dir_out;
                            T_light *= exp(-sigma_m * dt);
                            p_trans_dir *= exp(-sigma_m * dt);
                            p_trans_nee *= exp(-sigma_m * dt);
                            break;
                        }
                        iteration++;
                    }
                }
                if (!intersection_) {
                    break;
                } else {
                    PathVertex intersection = *intersection_;
                    if (intersection.material_id >= 0) {
                        T_light = make_zero_spectrum();
                        break;
                    }
                    shadow_bounces++;
                    if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                        T_light = make_zero_spectrum();
                        break;
                    }
                    if (intersection.material_id == -1) {
                        updateMedium(intersection, shadow_ray.dir, shadow_medium);
                    }
                    p = intersection.position;
                }
            }

            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);

            if (max(T_light) > 0) {
                T_light.x = T_light.x > 0 ? T_light.x : 0;
                T_light.y = T_light.y > 0 ? T_light.y : 0;
                T_light.z = T_light.z > 0 ? T_light.z : 0;
                Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p_org, scene);
                if (pdf_nee <= 0) {
                    break;
                }
                pdf_nee *= avg(p_trans_nee);
                Real G = abs(dot(dir_out, p_light.normal)) / distance_squared(ray.org, p_light.position);
                Spectrum Le_light = emission(light, -dir_out, Real(0), p_light, scene);
                Spectrum phase = eval(phase_f, dir_in, dir_out);
                Real pdf_phase = pdf_sample_phase(phase_f, dir_in, dir_out) * G * avg(p_trans_dir);
                Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
                Spectrum contribution = T_light * G * phase * Le_light / pdf_nee;
                radiance += current_path_throughput * sigma_s * contribution * w;
            }

            // Phase sampling
            Vector2 phase_uv = Vector2{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto next_dir_ = sample_phase_function(phase_f, -ray.dir, phase_uv);
            Vector3 next_dir;
            if (next_dir_) {
                next_dir = *next_dir_;
            }
            Spectrum phase = eval(phase_f, -ray.dir, next_dir);
            Real pdf_phase = pdf_sample_phase(phase_f, -ray.dir, next_dir);

            current_path_throughput *= (phase / pdf_phase) * sigma_s;

            ray.dir = next_dir;

            pos_cache = ray.org;
            pdf_sampling_cache = pdf_phase;
            pdf_trans_dir_cache = make_const_spectrum(1);
            pdf_trans_nee_cache = make_const_spectrum(1);
            never_scatter_or_bsdf = false;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(length(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    return radiance;
}
