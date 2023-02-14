#pragma once

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
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
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
        Real transmittance_pdf = exp(-sigma_t * t) * sigma_t;
        transmittance = exp(-sigma_t_v * t);

        // current point p
        Vector3 p = ray.org + t * ray.dir;

        // sample a point p' on light
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal p_light = sample_point_on_light(light, p, light_uv, shape_w, scene);

        Vector3 dir_out = normalize(p_light.position - p);
        Spectrum phase = eval(get_phase_function(scene.media[0]), -ray.dir, dir_out);  // only 1 volume
        Spectrum Le_light = emission(light,
                                     -dir_out,
                                     Real(0),
                                     p_light,
                                     scene);

        Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p, scene);

        Spectrum L_s1_estimate = phase * Le_light * exp(-sigma_t * distance(p, p_light.position)) * (abs(dot(dir_out, p_light.normal)) / distance_squared(p, p_light.position));

        return (transmittance / transmittance_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf);
    } else {
        // Real transmittance_pdf = exp(-sigma_t * t_hit) + sigma_t;
        Real transmittance_pdf = exp(-sigma_t * t_hit);
        transmittance = exp(-sigma_t_v * t_hit);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        return (transmittance / transmittance_pdf) * Le;
    }

    return make_zero_spectrum();
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
        Real transmittance_pdf = 1;

        if (current_medium >= 0) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t_v = sigma_a + sigma_s;
            Real sigma_t = sigma_a.x + sigma_s.x;  // the volume is monochromatic
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            if (t < t_hit) {
                transmittance_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t_v * t);
                scatter = true;
                ray.org = ray.org + t * ray.dir;
            } else {
                transmittance_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t_v * t_hit);
            }
        }

        current_path_throughput *= transmittance / transmittance_pdf;

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
        Real transmittance_pdf = 1;

        if (current_medium >= 0) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t_v = sigma_a + sigma_s;
            Real sigma_t = sigma_a.x + sigma_s.x;  // the volume is monochromatic
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            if (t < t_hit) {
                transmittance_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t_v * t);
                scatter = true;
                never_scatter = false;
                ray.org = ray.org + t * ray.dir;
            } else {
                transmittance_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t_v * t_hit);
                ray.org = ray.org + t_hit * ray.dir;
            }
            pdf_multi_trans_cache *= transmittance_pdf;
        }

        current_path_throughput *= transmittance / transmittance_pdf;

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
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            int light_id = sample_light(scene, light_w);
            const Light &light = scene.lights[light_id];
            PointAndNormal p_light = sample_point_on_light(light, ray.org, light_uv, shape_w, scene);

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
    Real pdf_sampling_phase_cache = 0;
    Real pdf_sampling_bsdf_cache = 0;
    Vector3 pos_cache = ray.org;  // for NEE
    Spectrum phase_cache = make_zero_spectrum();
    Real pdf_multi_trans_cache = 1;
    bool never_scatter_or_bsdf = true;
    bool last_is_scatter = false;

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
        Real transmittance_pdf = 1;

        if (current_medium >= 0) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t_v = sigma_a + sigma_s;
            Real sigma_t = sigma_a.x + sigma_s.x;  // the volume is monochromatic
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            if (t < t_hit) {
                transmittance_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t_v * t);
                scatter = true;
                never_scatter_or_bsdf = false;
                ray.org = ray.org + t * ray.dir;
            } else {
                transmittance_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t_v * t_hit);
                ray.org = ray.org + t_hit * ray.dir;
            }
            pdf_multi_trans_cache *= transmittance_pdf;
        }

        current_path_throughput *= transmittance / transmittance_pdf;

        // hit object
        if (!scatter && vertex.shape_id >= 0) {
            if (never_scatter_or_bsdf) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    Spectrum Le = emission(vertex, -ray.dir, scene);
                    radiance += current_path_throughput * Le;
                }
            } else {
                // need to consider MIS
                if (last_is_scatter) {
                    if (is_light(scene.shapes[vertex.shape_id])) {
                        int id = get_area_light_id(scene.shapes[vertex.shape_id]);
                        Real pdf_nee = light_pmf(scene, id) * pdf_point_on_light(scene.lights[id],
                                                                                 PointAndNormal{vertex.position, vertex.geometric_normal},
                                                                                 pos_cache,
                                                                                 scene);
                        Vector3 dir_out = normalize(vertex.position - pos_cache);
                        Real G = abs(dot(dir_out, vertex.geometric_normal)) / distance_squared(pos_cache, vertex.position);
                        Real pdf_dir = pdf_sampling_phase_cache * pdf_multi_trans_cache * G;
                        Real w_phase = (pdf_dir * pdf_dir) / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                        radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w_phase;
                    }
                } else {
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
                    if (is_light(scene.shapes[vertex.shape_id])) {
                        int id = get_area_light_id(scene.shapes[vertex.shape_id]);
                        Real pdf_nee = light_pmf(scene, id) * pdf_point_on_light(scene.lights[id],
                                                                                 PointAndNormal{vertex.position, vertex.geometric_normal},
                                                                                 pos_cache,
                                                                                 scene);
                        Vector3 dir_out = normalize(vertex.position - pos_cache);
                        Real G = abs(dot(dir_out, vertex.geometric_normal)) / distance_squared(pos_cache, vertex.position);
                        Real pdf_dir = pdf_sampling_bsdf_cache * pdf_multi_trans_cache * G;
                        Real w_bsdf = (pdf_dir * pdf_dir) / (pdf_dir * pdf_dir + pdf_nee * pdf_nee);
                        radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w_bsdf;
                    }
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
            // Volume Scatter
            // --------------------------------

            // NEE
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            int light_id = sample_light(scene, light_w);
            const Light &light = scene.lights[light_id];
            PointAndNormal p_light = sample_point_on_light(light, ray.org, light_uv, shape_w, scene);

            Vector3 p = ray.org;
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

            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            PhaseFunction phase_f = get_phase_function(scene.media[current_medium]);

            if (T_light > 0) {
                Vector3 dir_out = normalize(p_light.position - ray.org);
                Real G = abs(dot(dir_out, p_light.normal)) / distance_squared(ray.org, p_light.position);
                Spectrum Le_light = emission(light, -dir_out, Real(0), p_light, scene);
                Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, ray.org, scene);
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
            last_is_scatter = true;
        } else {
            // Hit Surface
            // ---------------------------------
            const Material &mat = scene.materials[vertex.material_id];
            const Vector3 &dir_in = -ray.dir;

            // NEE
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            int light_id = sample_light(scene, light_w);
            const Light &light = scene.lights[light_id];
            PointAndNormal p_light = sample_point_on_light(light, ray.org, light_uv, shape_w, scene);

            Vector3 p_org = vertex.position;
            Vector3 p = p_org;
            Vector3 dir_nee = normalize(p_light.position - p_org);
            Real T_light = 1;
            int shadow_medium = current_medium;
            int shadow_bounces = 0;
            Real p_trans_dir = 1;
            while (true) {
                Ray shadow_ray = Ray{p, dir_nee, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_light.position, p)};
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

            if (T_light > 0) {
                Vector3 dir_out = normalize(p_light.position - p_org);
                Real G = abs(dot(dir_out, p_light.normal)) / distance_squared(p_org, p_light.position);
                Spectrum Le_light = emission(light, -dir_out, Real(0), p_light, scene);
                Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_light, p_org, scene);
                Spectrum f = eval(mat, dir_in, dir_nee, vertex, scene.texture_pool);
                Real pdf_bsdf = G * pdf_sample_bsdf(mat, dir_in, dir_nee, vertex, scene.texture_pool);
                Spectrum contribution = T_light * G * f * Le_light / pdf_nee;
                Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
                radiance += current_path_throughput * contribution * w;
            }

            // BSDF sampling
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ = sample_bsdf(mat, dir_in, vertex, scene.texture_pool, bsdf_rnd_param_uv, bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }
            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
            Vector3 next_dir = bsdf_sample.dir_out;
            Ray bsdf_ray{vertex.position, next_dir, get_intersection_epsilon(scene), infinity<Real>()};
            std::optional<PathVertex> bsdf_intersection = intersect(scene, bsdf_ray);
            Real G = 1;
            if (bsdf_intersection) {
                G = fabs(dot(next_dir, bsdf_intersection->geometric_normal)) / distance_squared(bsdf_intersection->position, vertex.position);
            }
            Spectrum bsdf = eval(mat, dir_in, next_dir, vertex, scene.texture_pool);
            Real pdf_bsdf = pdf_sample_bsdf(mat, dir_in, next_dir, vertex, scene.texture_pool);
            if (pdf_bsdf <= 0) {
                // Numerical issue -- we generated some invalid rays.
                break;
            }
            pdf_sampling_bsdf_cache = pdf_bsdf;
            // pdf_bsdf *= G;
            current_path_throughput *= bsdf / pdf_bsdf;

            pos_cache = ray.org;
            pdf_multi_trans_cache = 1;
            ray.org = vertex.position;
            ray.dir = next_dir;
            last_is_scatter = false;
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
    return make_zero_spectrum();
}
