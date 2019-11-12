#include <jni.h>
#include <string>
#include <iostream>
#include <sstream>
#include <android/log.h>

#include <boost/filesystem/operations.hpp>   // for create_directories, exists
#include <boost/filesystem/path.hpp>         // for path, operator<<
#include <boost/filesystem/path_traits.hpp>  // for filesystem
#include <stddef.h>                          // for size_t
#include <sys/stat.h>                        // for stat
#include <volk/volk_prefs.h>                 // for volk_get_config_path
#include <map>                               // for map, map<>::iterator
#include <utility>                           // for pair

#include <volk/volk.h>                              // for volk_func_desc_t
#include <volk/volk_malloc.h>                       // for volk_free, volk_m...

#include <assert.h>                                 // for assert
#include <stdint.h>                                 // for uint16_t, uint64_t
#include <sys/time.h>                               // for CLOCKS_PER_SEC
#include <sys/types.h>                              // for int16_t, int32_t
#include <chrono>
#include <cmath>                                    // for sqrt, fabs, abs
#include <cstring>                                  // for memcpy, memset
#include <ctime>                                    // for clock
#include <fstream>                                  // for operator<<, basic...
#include <iostream>                                 // for cout, cerr
#include <limits>                                   // for numeric_limits
#include <map>                                      // for map, map<>::mappe...
#include <random>
#include <vector>                                   // for vector, _Bit_refe...



struct volk_type_t {
    bool is_float;
    bool is_scalar;
    bool is_signed;
    bool is_complex;
    int size;
    std::string str;
};

class volk_test_time_t {
public:
    std::string name;
    double time;
    std::string units;
    bool pass;
};

class volk_test_results_t {
public:
    std::string name;
    std::string config_name;
    unsigned int vlen;
    unsigned int iter;
    std::map<std::string, volk_test_time_t> results;
    std::string best_arch_a;
    std::string best_arch_u;
};

class volk_test_params_t {
private:
    float _tol;
    lv_32fc_t _scalar;
    unsigned int _vlen;
    unsigned int _iter;
    bool _benchmark_mode;
    bool _absolute_mode;
    std::string _kernel_regex;
public:
    // ctor
    volk_test_params_t(float tol, lv_32fc_t scalar, unsigned int vlen, unsigned int iter,
                       bool benchmark_mode, std::string kernel_regex) :
            _tol(tol), _scalar(scalar), _vlen(vlen), _iter(iter),
            _benchmark_mode(benchmark_mode), _absolute_mode(false), _kernel_regex(kernel_regex) {};
    // setters
    void set_tol(float tol) {_tol=tol;};
    void set_scalar(lv_32fc_t scalar) {_scalar=scalar;};
    void set_vlen(unsigned int vlen) {_vlen=vlen;};
    void set_iter(unsigned int iter) {_iter=iter;};
    void set_benchmark(bool benchmark) {_benchmark_mode=benchmark;};
    void set_regex(std::string regex) {_kernel_regex=regex;};
    // getters
    float tol() {return _tol;};
    lv_32fc_t scalar() {return _scalar;};
    unsigned int vlen() {return _vlen;};
    unsigned int iter() {return _iter;};
    bool benchmark_mode() {return _benchmark_mode;};
    bool absolute_mode() {return _absolute_mode;};
    std::string kernel_regex() {return _kernel_regex;};
    volk_test_params_t make_absolute(float tol) {
        volk_test_params_t t(*this);
        t._tol = tol;
        t._absolute_mode = true;
        return t;
    }
    volk_test_params_t make_tol(float tol) {
        volk_test_params_t t(*this);
        t._tol = tol;
        return t;
    }
};

class volk_test_case_t {
private:
    volk_func_desc_t _desc;
    void(*_kernel_ptr)();
    std::string _name;
    volk_test_params_t _test_parameters;
    std::string _puppet_master_name;
public:
    volk_func_desc_t desc() {return _desc;};
    void (*kernel_ptr()) () {return _kernel_ptr;};
    std::string name() {return _name;};
    std::string puppet_master_name() {return _puppet_master_name;};
    volk_test_params_t test_parameters() {return _test_parameters;};
    // normal ctor
    volk_test_case_t(volk_func_desc_t desc, void(*kernel_ptr)(), std::string name,
                     volk_test_params_t test_parameters) :
            _desc(desc), _kernel_ptr(kernel_ptr), _name(name), _test_parameters(test_parameters),
            _puppet_master_name("NULL")
    {};
    // ctor for puppets
    volk_test_case_t(volk_func_desc_t desc, void(*kernel_ptr)(), std::string name,
                     std::string puppet_master_name, volk_test_params_t test_parameters) :
            _desc(desc), _kernel_ptr(kernel_ptr), _name(name), _test_parameters(test_parameters),
            _puppet_master_name(puppet_master_name)
    {};
};

/************************************************
 * VOLK QA functions                            *
 ************************************************/
volk_type_t volk_type_from_string(std::string);

float uniform(void);
void random_floats(float *buf, unsigned n);

bool run_volk_tests(
        volk_func_desc_t,
        void(*)(),
        std::string,
        volk_test_params_t,
        std::vector<volk_test_results_t> *results = NULL,
        std::string puppet_master_name = "NULL"
);

bool run_volk_tests(
        volk_func_desc_t,
        void(*)(),
        std::string,
        float,
        lv_32fc_t,
        unsigned int,
        unsigned int,
        std::vector<volk_test_results_t> *results = NULL,
        std::string puppet_master_name = "NULL",
        bool absolute_mode = false,
        bool benchmark_mode = false
);


#define VOLK_RUN_TESTS(func, tol, scalar, len, iter) \
    BOOST_AUTO_TEST_CASE(func##_test) { \
        BOOST_CHECK_EQUAL(run_volk_tests( \
            func##_get_func_desc(), (void (*)())func##_manual, \
            std::string(#func), tol, scalar, len, iter, 0, "NULL"), \
          0); \
    }
#define VOLK_PROFILE(func, test_params, results) run_volk_tests(func##_get_func_desc(), (void (*)())func##_manual, std::string(#func), test_params, results, "NULL")
#define VOLK_PUPPET_PROFILE(func, puppet_master_func, test_params, results) run_volk_tests(func##_get_func_desc(), (void (*)())func##_manual, std::string(#func), test_params, results, std::string(#puppet_master_func))
typedef void (*volk_fn_1arg)(void *, unsigned int, const char*); //one input, operate in place
typedef void (*volk_fn_2arg)(void *, void *, unsigned int, const char*);
typedef void (*volk_fn_3arg)(void *, void *, void *, unsigned int, const char*);
typedef void (*volk_fn_4arg)(void *, void *, void *, void *, unsigned int, const char*);
typedef void (*volk_fn_1arg_s32f)(void *, float, unsigned int, const char*); //one input vector, one scalar float input
typedef void (*volk_fn_2arg_s32f)(void *, void *, float, unsigned int, const char*);
typedef void (*volk_fn_3arg_s32f)(void *, void *, void *, float, unsigned int, const char*);
typedef void (*volk_fn_1arg_s32fc)(void *, lv_32fc_t, unsigned int, const char*); //one input vector, one scalar float input
typedef void (*volk_fn_2arg_s32fc)(void *, void *, lv_32fc_t, unsigned int, const char*);
typedef void (*volk_fn_3arg_s32fc)(void *, void *, void *, lv_32fc_t, unsigned int, const char*);












/************************************************
 * VOLK QA functions                            *
 ************************************************/
volk_type_t volk_type_from_string(std::string);

float uniform(void);
void random_floats(float *buf, unsigned n);

#define VOLK_INIT_PUPP(func, puppet_master_func, test_params)\
    volk_test_case_t(func##_get_func_desc(), (void(*)())func##_manual, std::string(#func),\
    std::string(#puppet_master_func), test_params)

#define VOLK_INIT_TEST(func, test_params)\
    volk_test_case_t(func##_get_func_desc(), (void(*)())func##_manual, std::string(#func),\
    test_params)

#define QA(test) test_cases.push_back(test);
std::vector<volk_test_case_t> init_test_list(volk_test_params_t test_params)
{

    // Some kernels need a lower tolerance
    volk_test_params_t test_params_inacc = test_params.make_tol(1e-2);
    volk_test_params_t test_params_inacc_tenth = test_params.make_tol(1e-1);

    volk_test_params_t test_params_power(test_params);
    test_params_power.set_scalar(2.5);

    volk_test_params_t test_params_rotator(test_params);
    test_params_rotator.set_scalar(std::polar(1.0f, 0.1f));
    test_params_rotator.set_tol(1e-3);

    std::vector<volk_test_case_t> test_cases;
    QA(VOLK_INIT_PUPP(volk_64u_popcntpuppet_64u, volk_64u_popcnt,     test_params))
    QA(VOLK_INIT_PUPP(volk_64u_popcntpuppet_64u, volk_64u_popcnt,     test_params))
    QA(VOLK_INIT_PUPP(volk_64u_popcntpuppet_64u, volk_64u_popcnt,     test_params))
    QA(VOLK_INIT_PUPP(volk_16u_byteswappuppet_16u, volk_16u_byteswap, test_params))
    QA(VOLK_INIT_PUPP(volk_32u_byteswappuppet_32u, volk_32u_byteswap, test_params))
    QA(VOLK_INIT_PUPP(volk_32u_popcntpuppet_32u, volk_32u_popcnt_32u,  test_params))
    QA(VOLK_INIT_PUPP(volk_64u_byteswappuppet_64u, volk_64u_byteswap, test_params))
    QA(VOLK_INIT_PUPP(volk_32fc_s32fc_rotatorpuppet_32fc, volk_32fc_s32fc_x2_rotator_32fc, test_params_rotator))
    QA(VOLK_INIT_PUPP(volk_8u_conv_k7_r2puppet_8u, volk_8u_x4_conv_k7_r2_8u, test_params.make_tol(0)))
    QA(VOLK_INIT_PUPP(volk_32f_x2_fm_detectpuppet_32f, volk_32f_s32f_32f_fm_detect_32f, test_params))
    QA(VOLK_INIT_TEST(volk_16ic_s32f_deinterleave_real_32f,           test_params))
    QA(VOLK_INIT_TEST(volk_16ic_deinterleave_real_8i,                 test_params))
    QA(VOLK_INIT_TEST(volk_16ic_deinterleave_16i_x2,                  test_params))
    QA(VOLK_INIT_TEST(volk_16ic_s32f_deinterleave_32f_x2,             test_params))
    QA(VOLK_INIT_TEST(volk_16ic_deinterleave_real_16i,                test_params))
    QA(VOLK_INIT_TEST(volk_16ic_magnitude_16i,                        test_params))
    QA(VOLK_INIT_TEST(volk_16ic_s32f_magnitude_32f,                   test_params))
    QA(VOLK_INIT_TEST(volk_16ic_convert_32fc,                         test_params))
    QA(VOLK_INIT_TEST(volk_16ic_x2_multiply_16ic,                     test_params))
    QA(VOLK_INIT_TEST(volk_16ic_x2_dot_prod_16ic,                     test_params))
    QA(VOLK_INIT_TEST(volk_16i_s32f_convert_32f,                      test_params))
    QA(VOLK_INIT_TEST(volk_16i_convert_8i,                            test_params))
    QA(VOLK_INIT_TEST(volk_16i_32fc_dot_prod_32fc,                    test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_accumulator_s32f,                      test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_x2_add_32f,                            test_params))
    QA(VOLK_INIT_TEST(volk_32f_index_max_16u,                         test_params))
    QA(VOLK_INIT_TEST(volk_32f_index_max_32u,                         test_params))
    QA(VOLK_INIT_TEST(volk_32fc_32f_multiply_32fc,                    test_params))
    QA(VOLK_INIT_TEST(volk_32fc_32f_add_32fc,                         test_params))
    QA(VOLK_INIT_TEST(volk_32f_log2_32f,                              test_params.make_absolute(1e-5)))
    QA(VOLK_INIT_TEST(volk_32f_expfast_32f,                           test_params_inacc_tenth))
    QA(VOLK_INIT_TEST(volk_32f_x2_pow_32f,                            test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_sin_32f,                               test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_cos_32f,                               test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_tan_32f,                               test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_atan_32f,                              test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_asin_32f,                              test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_acos_32f,                              test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32fc_s32f_power_32fc,                      test_params_power))
    QA(VOLK_INIT_TEST(volk_32f_s32f_calc_spectral_noise_floor_32f,    test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32fc_s32f_atan2_32f,                       test_params))
    QA(VOLK_INIT_TEST(volk_32fc_x2_conjugate_dot_prod_32fc,           test_params_inacc_tenth))
    QA(VOLK_INIT_TEST(volk_32fc_deinterleave_32f_x2,                  test_params))
    QA(VOLK_INIT_TEST(volk_32fc_deinterleave_64f_x2,                  test_params))
    QA(VOLK_INIT_TEST(volk_32fc_s32f_deinterleave_real_16i,           test_params))
    QA(VOLK_INIT_TEST(volk_32fc_deinterleave_imag_32f,                test_params))
    QA(VOLK_INIT_TEST(volk_32fc_deinterleave_real_32f,                test_params))
    QA(VOLK_INIT_TEST(volk_32fc_deinterleave_real_64f,                test_params))
    QA(VOLK_INIT_TEST(volk_32fc_x2_dot_prod_32fc,                     test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32fc_32f_dot_prod_32fc,                    test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32fc_index_max_16u,                        test_params))
    QA(VOLK_INIT_TEST(volk_32fc_index_max_32u,                        test_params))
    QA(VOLK_INIT_TEST(volk_32fc_s32f_magnitude_16i,                   test_params))
    QA(VOLK_INIT_TEST(volk_32fc_magnitude_32f,                        test_params_inacc_tenth))
    QA(VOLK_INIT_TEST(volk_32fc_magnitude_squared_32f,                test_params))
    QA(VOLK_INIT_TEST(volk_32fc_x2_add_32fc,                          test_params))
    QA(VOLK_INIT_TEST(volk_32fc_x2_multiply_32fc,                     test_params))
    QA(VOLK_INIT_TEST(volk_32fc_x2_multiply_conjugate_32fc,           test_params))
    QA(VOLK_INIT_TEST(volk_32fc_x2_divide_32fc,                       test_params))
    QA(VOLK_INIT_TEST(volk_32fc_conjugate_32fc,                       test_params))
    QA(VOLK_INIT_TEST(volk_32f_s32f_convert_16i,                      test_params))
    QA(VOLK_INIT_TEST(volk_32f_s32f_convert_32i,                      test_params))
    QA(VOLK_INIT_TEST(volk_32f_convert_64f,                           test_params))
    QA(VOLK_INIT_TEST(volk_32f_s32f_convert_8i,                       test_params))
    QA(VOLK_INIT_TEST(volk_32fc_convert_16ic,                         test_params))
    QA(VOLK_INIT_TEST(volk_32fc_s32f_power_spectrum_32f,              test_params))
    QA(VOLK_INIT_TEST(volk_32fc_x2_square_dist_32f,                   test_params))
    QA(VOLK_INIT_TEST(volk_32fc_x2_s32f_square_dist_scalar_mult_32f,  test_params))
    QA(VOLK_INIT_TEST(volk_32f_x2_divide_32f,                         test_params))
    QA(VOLK_INIT_TEST(volk_32f_x2_dot_prod_32f,                       test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_x2_s32f_interleave_16ic,               test_params))
    QA(VOLK_INIT_TEST(volk_32f_x2_interleave_32fc,                    test_params))
    QA(VOLK_INIT_TEST(volk_32f_x2_max_32f,                            test_params))
    QA(VOLK_INIT_TEST(volk_32f_x2_min_32f,                            test_params))
    QA(VOLK_INIT_TEST(volk_32f_x2_multiply_32f,                       test_params))
    QA(VOLK_INIT_TEST(volk_32f_64f_multiply_64f,                      test_params))
    QA(VOLK_INIT_TEST(volk_32f_64f_add_64f,                           test_params))
    QA(VOLK_INIT_TEST(volk_32f_s32f_normalize,                        test_params))
    QA(VOLK_INIT_TEST(volk_32f_s32f_power_32f,                        test_params))
    QA(VOLK_INIT_TEST(volk_32f_sqrt_32f,                              test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_s32f_stddev_32f,                       test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_stddev_and_mean_32f_x2,                test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_x2_subtract_32f,                       test_params))
    QA(VOLK_INIT_TEST(volk_32f_x3_sum_of_poly_32f,                    test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32i_x2_and_32i,                            test_params))
    QA(VOLK_INIT_TEST(volk_32i_s32f_convert_32f,                      test_params))
    QA(VOLK_INIT_TEST(volk_32i_x2_or_32i,                             test_params))
    QA(VOLK_INIT_TEST(volk_32f_x2_dot_prod_16i,                       test_params))
    QA(VOLK_INIT_TEST(volk_64f_convert_32f,                           test_params))
    QA(VOLK_INIT_TEST(volk_64f_x2_max_64f,                            test_params))
    QA(VOLK_INIT_TEST(volk_64f_x2_min_64f,                            test_params))
    QA(VOLK_INIT_TEST(volk_64f_x2_multiply_64f,                       test_params))
    QA(VOLK_INIT_TEST(volk_64f_x2_add_64f,                            test_params))
    QA(VOLK_INIT_TEST(volk_8ic_deinterleave_16i_x2,                   test_params))
    QA(VOLK_INIT_TEST(volk_8ic_s32f_deinterleave_32f_x2,              test_params))
    QA(VOLK_INIT_TEST(volk_8ic_deinterleave_real_16i,                 test_params))
    QA(VOLK_INIT_TEST(volk_8ic_s32f_deinterleave_real_32f,            test_params))
    QA(VOLK_INIT_TEST(volk_8ic_deinterleave_real_8i,                  test_params))
    QA(VOLK_INIT_TEST(volk_8ic_x2_multiply_conjugate_16ic,            test_params))
    QA(VOLK_INIT_TEST(volk_8ic_x2_s32f_multiply_conjugate_32fc,       test_params))
    QA(VOLK_INIT_TEST(volk_8i_convert_16i,                            test_params))
    QA(VOLK_INIT_TEST(volk_8i_s32f_convert_32f,                       test_params))
    QA(VOLK_INIT_TEST(volk_32fc_s32fc_multiply_32fc,                  test_params))
    QA(VOLK_INIT_TEST(volk_32f_s32f_multiply_32f,                     test_params))
    QA(VOLK_INIT_TEST(volk_32f_binary_slicer_32i,                     test_params))
    QA(VOLK_INIT_TEST(volk_32f_binary_slicer_8i,                      test_params))
    QA(VOLK_INIT_TEST(volk_32u_reverse_32u,                            test_params))
    QA(VOLK_INIT_TEST(volk_32f_tanh_32f,                              test_params_inacc))
    QA(VOLK_INIT_TEST(volk_32f_s32f_mod_rangepuppet_32f,              test_params))
    QA(VOLK_INIT_PUPP(volk_8u_x3_encodepolarpuppet_8u, volk_8u_x3_encodepolar_8u_x2, test_params))
    QA(VOLK_INIT_PUPP(volk_32f_8u_polarbutterflypuppet_32f, volk_32f_8u_polarbutterfly_32f, test_params))
    // no one uses these, so don't test them
    //VOLK_PROFILE(volk_16i_x5_add_quad_16i_x4, 1e-4, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_16i_branch_4_state_8, 1e-4, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_16i_max_star_16i, 0, 0, 204602, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_16i_max_star_horizontal_16i, 0, 0, 204602, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_16i_permute_and_scalar_add, 1e-4, 0, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_16i_x4_quad_max_star_16i, 1e-4, 0, 2046, 10000, &results, benchmark_mode, kernel_regex);
    // we need a puppet for this one
    //(VOLK_INIT_TEST(volk_32fc_s32f_x2_power_spectral_density_32f,   test_params))


    return test_cases;
}









template <typename T>
void random_floats(void *buf, unsigned int n, std::default_random_engine& rnd_engine)
{
    T *array = static_cast<T*>(buf);
    std::uniform_real_distribution<T> uniform_dist(T(-1), T(1));
    for(unsigned int i = 0; i < n; i++) {
        array[i] = uniform_dist(rnd_engine);
    }
}

void load_random_data(void *data, volk_type_t type, unsigned int n) {
    std::random_device rnd_device;
    std::default_random_engine rnd_engine(rnd_device());
    if(type.is_complex) n *= 2;
    if(type.is_float) {
        if(type.size == 8) {
            random_floats<double>(data, n, rnd_engine);
        } else {
            random_floats<float> (data, n, rnd_engine);
        }
    } else {
        float int_max = float(uint64_t(2) << (type.size*8));
        if(type.is_signed) int_max /= 2.0;
        std::uniform_real_distribution<float> uniform_dist(-int_max, int_max);
        for(unsigned int i=0; i<n; i++) {
            float scaled_rand = uniform_dist(rnd_engine);
            //man i really don't know how to do this in a more clever way, you have to cast down at some point
            switch(type.size) {
                case 8:
                    if(type.is_signed) ((int64_t *)data)[i] = (int64_t) scaled_rand;
                    else ((uint64_t *)data)[i] = (uint64_t) scaled_rand;
                    break;
                case 4:
                    if(type.is_signed) ((int32_t *)data)[i] = (int32_t) scaled_rand;
                    else ((uint32_t *)data)[i] = (uint32_t) scaled_rand;
                    break;
                case 2:
                    if(type.is_signed) ((int16_t *)data)[i] = (int16_t)((int16_t) scaled_rand % 8);
                    else ((uint16_t *)data)[i] = (uint16_t) ((int16_t) scaled_rand % 8);
                    break;
                case 1:
                    if(type.is_signed) ((int8_t *)data)[i] = (int8_t) scaled_rand;
                    else ((uint8_t *)data)[i] = (uint8_t) scaled_rand;
                    break;
                default:
                    throw "load_random_data: no support for data size > 8 or < 1"; //no shenanigans here
            }
        }
    }
}

static std::vector<std::string> get_arch_list(volk_func_desc_t desc) {
    std::vector<std::string> archlist;

    for(size_t i = 0; i < desc.n_impls; i++) {
        archlist.push_back(std::string(desc.impl_names[i]));
    }

    return archlist;
}

template <typename T>
T volk_lexical_cast(const std::string& str)
{
    for (unsigned int c_index = 0; c_index < str.size(); ++c_index) {
        if (str.at(c_index) < '0' || str.at(c_index) > '9') {
            throw "not all numbers!";
        }
    }
    T var;
    std::istringstream iss;
    iss.str(str);
    iss >> var;
    // deal with any error bits that may have been set on the stream
    return var;
}

volk_type_t volk_type_from_string(std::string name) {
    volk_type_t type;
    type.is_float = false;
    type.is_scalar = false;
    type.is_complex = false;
    type.is_signed = false;
    type.size = 0;
    type.str = name;

    if(name.size() < 2) {
        throw std::string("name too short to be a datatype");
    }

    //is it a scalar?
    if(name[0] == 's') {
        type.is_scalar = true;
        name = name.substr(1, name.size()-1);
    }

    //get the data size
    size_t last_size_pos = name.find_last_of("0123456789");
    if(last_size_pos == std::string::npos) {
        throw std::string("no size spec in type ").append(name);
    }
    //will throw if malformed
    int size = volk_lexical_cast<int>(name.substr(0, last_size_pos+1));

    assert(((size % 8) == 0) && (size <= 64) && (size != 0));
    type.size = size/8; //in bytes

    for(size_t i=last_size_pos+1; i < name.size(); i++) {
        switch (name[i]) {
            case 'f':
                type.is_float = true;
                break;
            case 'i':
                type.is_signed = true;
                break;
            case 'c':
                type.is_complex = true;
                break;
            case 'u':
                type.is_signed = false;
                break;
            default:
                throw;
        }
    }

    return type;
}

std::vector<std::string> split_signature(const std::string &protokernel_signature) {
    std::vector<std::string> signature_tokens;
    std::string token;
    for (unsigned int loc = 0; loc < protokernel_signature.size(); ++loc) {
        if (protokernel_signature.at(loc) == '_') {
            // this is a break
            signature_tokens.push_back(token);
            token = "";
        } else {
            token.push_back(protokernel_signature.at(loc));
        }
    }
    // Get the last one to the end of the string
    signature_tokens.push_back(token);
    return signature_tokens;
}

static void get_signatures_from_name(std::vector<volk_type_t> &inputsig,
                                     std::vector<volk_type_t> &outputsig,
                                     std::string name) {

    std::vector<std::string> toked = split_signature(name);

    assert(toked[0] == "volk");
    toked.erase(toked.begin());

    //ok. we're assuming a string in the form
    //(sig)_(multiplier-opt)_..._(name)_(sig)_(multiplier-opt)_..._(alignment)

    enum { SIDE_INPUT, SIDE_NAME, SIDE_OUTPUT } side = SIDE_INPUT;
    std::string fn_name;
    volk_type_t type;
    for (unsigned int token_index = 0; token_index < toked.size(); ++token_index) {
        std::string token = toked[token_index];
        try {
            type = volk_type_from_string(token);
            if(side == SIDE_NAME) side = SIDE_OUTPUT; //if this is the first one after the name...

            if(side == SIDE_INPUT) inputsig.push_back(type);
            else outputsig.push_back(type);
        } catch (...){
            if(token[0] == 'x' && (token.size() > 1) && (token[1] > '0' || token[1] < '9')) { //it's a multiplier
                if(side == SIDE_INPUT) assert(inputsig.size() > 0);
                else assert(outputsig.size() > 0);
                int multiplier = volk_lexical_cast<int>(token.substr(1, token.size()-1)); //will throw if invalid
                for(int i=1; i<multiplier; i++) {
                    if(side == SIDE_INPUT) inputsig.push_back(inputsig.back());
                    else outputsig.push_back(outputsig.back());
                }
            }
            else if(side == SIDE_INPUT) { //it's the function name, at least it better be
                side = SIDE_NAME;
                fn_name.append("_");
                fn_name.append(token);
            }
            else if(side == SIDE_OUTPUT) {
                if(token != toked.back()) throw; //the last token in the name is the alignment
            }
        }
    }
    //we don't need an output signature (some fn's operate on the input data, "in place"), but we do need at least one input!
    assert(inputsig.size() != 0);

}

inline void run_cast_test1(volk_fn_1arg func, std::vector<void *> &buffs, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], vlen, arch.c_str());
}

inline void run_cast_test2(volk_fn_2arg func, std::vector<void *> &buffs, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], buffs[1], vlen, arch.c_str());
}

inline void run_cast_test3(volk_fn_3arg func, std::vector<void *> &buffs, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], buffs[1], buffs[2], vlen, arch.c_str());
}

inline void run_cast_test4(volk_fn_4arg func, std::vector<void *> &buffs, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], buffs[1], buffs[2], buffs[3], vlen, arch.c_str());
}

inline void run_cast_test1_s32f(volk_fn_1arg_s32f func, std::vector<void *> &buffs, float scalar, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], scalar, vlen, arch.c_str());
}

inline void run_cast_test2_s32f(volk_fn_2arg_s32f func, std::vector<void *> &buffs, float scalar, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], buffs[1], scalar, vlen, arch.c_str());
}

inline void run_cast_test3_s32f(volk_fn_3arg_s32f func, std::vector<void *> &buffs, float scalar, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], buffs[1], buffs[2], scalar, vlen, arch.c_str());
}

inline void run_cast_test1_s32fc(volk_fn_1arg_s32fc func, std::vector<void *> &buffs, lv_32fc_t scalar, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], scalar, vlen, arch.c_str());
}

inline void run_cast_test2_s32fc(volk_fn_2arg_s32fc func, std::vector<void *> &buffs, lv_32fc_t scalar, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], buffs[1], scalar, vlen, arch.c_str());
}

inline void run_cast_test3_s32fc(volk_fn_3arg_s32fc func, std::vector<void *> &buffs, lv_32fc_t scalar, unsigned int vlen, unsigned int iter, std::string arch) {
    while(iter--) func(buffs[0], buffs[1], buffs[2], scalar, vlen, arch.c_str());
}

template <class t>
bool fcompare(t *in1, t *in2, unsigned int vlen, float tol, bool absolute_mode) {
    bool fail = false;
    int print_max_errs = 10;
    for(unsigned int i=0; i<vlen; i++) {
        if (absolute_mode) {
            if (fabs(((t *)(in1))[i] - ((t *)(in2))[i]) > tol) {
                fail=true;
                if(print_max_errs-- > 0) {
                    std::cout << "offset " << i << " in1: " << t(((t *)(in1))[i]) << " in2: " << t(((t *)(in2))[i]);
                    std::cout << " tolerance was: " << tol << std::endl;
                }
            }
        } else {
            // for very small numbers we'll see round off errors due to limited
            // precision. So a special test case...
            if(fabs(((t *)(in1))[i]) < 1e-30) {
                if( fabs( ((t *)(in2))[i] ) > tol )
                {
                    fail=true;
                    if(print_max_errs-- > 0) {
                        std::cout << "offset " << i << " in1: " << t(((t *)(in1))[i]) << " in2: " << t(((t *)(in2))[i]);
                        std::cout << " tolerance was: " << tol << std::endl;
                    }
                }
            }
                // the primary test is the percent different greater than given tol
            else if(fabs(((t *)(in1))[i] - ((t *)(in2))[i])/fabs(((t *)in1)[i]) > tol) {
                fail=true;
                if(print_max_errs-- > 0) {
                    std::cout << "offset " << i << " in1: " << t(((t *)(in1))[i]) << " in2: " << t(((t *)(in2))[i]);
                    std::cout << " tolerance was: " << tol << std::endl;
                }
            }
        }
    }

    return fail;
}

template <class t>
bool ccompare(t *in1, t *in2, unsigned int vlen, float tol, bool absolute_mode) {
    if (absolute_mode) {
        std::cout << "ccompare does not support absolute mode" << std::endl;
        return true;
    }
    bool fail = false;
    int print_max_errs = 10;
    for(unsigned int i=0; i<2*vlen; i+=2) {
        if (std::isnan(in1[i]) || std::isnan(in1[i+1]) || std::isnan(in2[i]) || std::isnan(in2[i+1])
            || std::isinf(in1[i]) || std::isinf(in1[i+1]) || std::isinf(in2[i]) || std::isinf(in2[i+1])) {
            fail=true;
            if(print_max_errs-- > 0) {
                std::cout << "offset " << i/2 << " in1: " << in1[i] << " + " << in1[i+1] << "j  in2: " << in2[i] << " + " << in2[i+1] << "j";
                std::cout << " tolerance was: " << tol << std::endl;
            }
        }
        t diff[2] = { in1[i] - in2[i], in1[i+1] - in2[i+1] };
        t err  = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
        t norm = std::sqrt(in1[i] * in1[i] + in1[i+1] * in1[i+1]);

        // for very small numbers we'll see round off errors due to limited
        // precision. So a special test case...
        if (norm < 1e-30) {
            if (err > tol)
            {
                fail=true;
                if(print_max_errs-- > 0) {
                    std::cout << "offset " << i/2 << " in1: " << in1[i] << " + " << in1[i+1] << "j  in2: " << in2[i] << " + " << in2[i+1] << "j";
                    std::cout << " tolerance was: " << tol << std::endl;
                }
            }
        }
            // the primary test is the percent different greater than given tol
        else if((err / norm) > tol) {
            fail=true;
            if(print_max_errs-- > 0) {
                std::cout << "offset " << i/2 << " in1: " << in1[i] << " + " << in1[i+1] << "j  in2: " << in2[i] << " + " << in2[i+1] << "j";
                std::cout << " tolerance was: " << tol << std::endl;
            }
        }
    }

    return fail;
}

template <class t>
bool icompare(t *in1, t *in2, unsigned int vlen, unsigned int tol, bool absolute_mode) {
    if (absolute_mode) {
        std::cout << "icompare does not support absolute mode" << std::endl;
        return true;
    }
    bool fail = false;
    int print_max_errs = 10;
    for(unsigned int i=0; i<vlen; i++) {
        if(((unsigned int)abs(int(((t *)(in1))[i]) - int(((t *)(in2))[i]))) > tol) {
            fail=true;
            if(print_max_errs-- > 0) {
                std::cout << "offset " << i << " in1: " << static_cast<int>(t(((t *)(in1))[i])) << " in2: " << static_cast<int>(t(((t *)(in2))[i]));
                std::cout << " tolerance was: " << tol << std::endl;
            }
        }
    }

    return fail;
}

class volk_qa_aligned_mem_pool{
public:
    void *get_new(size_t size){
        size_t alignment = volk_get_alignment();
        void* ptr = volk_malloc(size, alignment);
        memset(ptr, 0x00, size);
        _mems.push_back(ptr);
        return ptr;
    }
    ~volk_qa_aligned_mem_pool() {
        for(unsigned int ii = 0; ii < _mems.size(); ++ii) {
            volk_free(_mems[ii]);
        }
    }
private: std::vector<void * > _mems;
};

bool run_volk_tests(volk_func_desc_t desc,
                    void (*manual_func)(),
                    std::string name,
                    volk_test_params_t test_params,
                    std::vector<volk_test_results_t> *results,
                    std::string puppet_master_name
)
{
    return run_volk_tests(desc, manual_func, name, test_params.tol(), test_params.scalar(),
                          test_params.vlen(), test_params.iter(), results, puppet_master_name,
                          test_params.absolute_mode(), test_params.benchmark_mode());
}

bool run_volk_tests(volk_func_desc_t desc,
                    void (*manual_func)(),
                    std::string name,
                    float tol,
                    lv_32fc_t scalar,
                    unsigned int vlen,
                    unsigned int iter,
                    std::vector<volk_test_results_t> *results,
                    std::string puppet_master_name,
                    bool absolute_mode,
                    bool benchmark_mode
) {
    // Initialize this entry in results vector
    results->push_back(volk_test_results_t());
    results->back().name = name;
    results->back().vlen = vlen;
    results->back().iter = iter;
    std::cout << "RUN_VOLK_TESTS: " << name << "(" << vlen << "," << iter << ")" << std::endl;

    // vlen_twiddle will increase vlen for malloc and data generation
    // but kernels will still be called with the user provided vlen.
    // This is useful for causing errors in kernels that do bad reads
    const unsigned int vlen_twiddle = 5;
    vlen = vlen + vlen_twiddle;

    const float tol_f = tol;
    const unsigned int tol_i = static_cast<const unsigned int>(tol);

    //first let's get a list of available architectures for the test
    std::vector<std::string> arch_list = get_arch_list(desc);

    if((!benchmark_mode) && (arch_list.size() < 2)) {
        std::cout << "no architectures to test" << std::endl;
        return false;
    }

    //something that can hang onto memory and cleanup when this function exits
    volk_qa_aligned_mem_pool mem_pool;

    //now we have to get a function signature by parsing the name
    std::vector<volk_type_t> inputsig, outputsig;
    try {
        get_signatures_from_name(inputsig, outputsig, name);
    }
    catch (std::exception &error) {
        std::cerr << "Error: unable to get function signature from kernel name" << std::endl;
        std::cerr << "  - " << name << std::endl;
        return false;
    }

    //pull the input scalars into their own vector
    std::vector<volk_type_t> inputsc;
    for(size_t i=0; i<inputsig.size(); i++) {
        if(inputsig[i].is_scalar) {
            inputsc.push_back(inputsig[i]);
            inputsig.erase(inputsig.begin() + i);
            i -= 1;
        }
    }
    std::vector<void *> inbuffs;
    for (unsigned int inputsig_index = 0; inputsig_index < inputsig.size(); ++ inputsig_index) {
        volk_type_t sig = inputsig[inputsig_index];
        if(!sig.is_scalar) //we don't make buffers for scalars
            inbuffs.push_back(mem_pool.get_new(vlen*sig.size*(sig.is_complex ? 2 : 1)));
    }
    for(size_t i=0; i<inbuffs.size(); i++) {
        load_random_data(inbuffs[i], inputsig[i], vlen);
    }

    //ok let's make a vector of vector of void buffers, which holds the input/output vectors for each arch
    std::vector<std::vector<void *> > test_data;
    for(size_t i=0; i<arch_list.size(); i++) {
        std::vector<void *> arch_buffs;
        for(size_t j=0; j<outputsig.size(); j++) {
            arch_buffs.push_back(mem_pool.get_new(vlen*outputsig[j].size*(outputsig[j].is_complex ? 2 : 1)));
        }
        for(size_t j=0; j<inputsig.size(); j++) {
            void *arch_inbuff = mem_pool.get_new(vlen*inputsig[j].size*(inputsig[j].is_complex ? 2 : 1));
            memcpy(arch_inbuff, inbuffs[j], vlen * inputsig[j].size * (inputsig[j].is_complex ? 2 : 1));
            arch_buffs.push_back(arch_inbuff);
        }
        test_data.push_back(arch_buffs);
    }

    std::vector<volk_type_t> both_sigs;
    both_sigs.insert(both_sigs.end(), outputsig.begin(), outputsig.end());
    both_sigs.insert(both_sigs.end(), inputsig.begin(), inputsig.end());

    //now run the test
    vlen = vlen - vlen_twiddle;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::vector<double> profile_times;
    for(size_t i = 0; i < arch_list.size(); i++) {
        start = std::chrono::system_clock::now();

        switch(both_sigs.size()) {
            case 1:
                if(inputsc.size() == 0) {
                    run_cast_test1((volk_fn_1arg)(manual_func), test_data[i], vlen, iter, arch_list[i]);
                } else if(inputsc.size() == 1 && inputsc[0].is_float) {
                    if(inputsc[0].is_complex) {
                        run_cast_test1_s32fc((volk_fn_1arg_s32fc)(manual_func), test_data[i], scalar, vlen, iter, arch_list[i]);
                    } else {
                        run_cast_test1_s32f((volk_fn_1arg_s32f)(manual_func), test_data[i], scalar.real(), vlen, iter, arch_list[i]);
                    }
                } else throw "unsupported 1 arg function >1 scalars";
                break;
            case 2:
                if(inputsc.size() == 0) {
                    run_cast_test2((volk_fn_2arg)(manual_func), test_data[i], vlen, iter, arch_list[i]);
                } else if(inputsc.size() == 1 && inputsc[0].is_float) {
                    if(inputsc[0].is_complex) {
                        run_cast_test2_s32fc((volk_fn_2arg_s32fc)(manual_func), test_data[i], scalar, vlen, iter, arch_list[i]);
                    } else {
                        run_cast_test2_s32f((volk_fn_2arg_s32f)(manual_func), test_data[i], scalar.real(), vlen, iter, arch_list[i]);
                    }
                } else throw "unsupported 2 arg function >1 scalars";
                break;
            case 3:
                if(inputsc.size() == 0) {
                    run_cast_test3((volk_fn_3arg)(manual_func), test_data[i], vlen, iter, arch_list[i]);
                } else if(inputsc.size() == 1 && inputsc[0].is_float) {
                    if(inputsc[0].is_complex) {
                        run_cast_test3_s32fc((volk_fn_3arg_s32fc)(manual_func), test_data[i], scalar, vlen, iter, arch_list[i]);
                    } else {
                        run_cast_test3_s32f((volk_fn_3arg_s32f)(manual_func), test_data[i], scalar.real(), vlen, iter, arch_list[i]);
                    }
                } else throw "unsupported 3 arg function >1 scalars";
                break;
            case 4:
                run_cast_test4((volk_fn_4arg)(manual_func), test_data[i], vlen, iter, arch_list[i]);
                break;
            default:
                throw "no function handler for this signature";
                break;
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        double arch_time = 1000.0 * elapsed_seconds.count();
        std::cout << arch_list[i] << " completed in " << arch_time << " ms" << std::endl;
        volk_test_time_t result;
        result.name = arch_list[i];
        result.time = arch_time;
        result.units = "ms";
        result.pass = true;
        results->back().results[result.name] = result;

        profile_times.push_back(arch_time);
    }

    //and now compare each output to the generic output
    //first we have to know which output is the generic one, they aren't in order...
    size_t generic_offset=0;
    for(size_t i=0; i<arch_list.size(); i++) {
        if (arch_list[i] == "generic") {
            generic_offset = i;
        }
    }

    // Just in case a kernel wrote to OOB memory, use the twiddled vlen
    vlen = vlen + vlen_twiddle;
    bool fail;
    bool fail_global = false;
    std::vector<bool> arch_results;
    for(size_t i=0; i<arch_list.size(); i++) {
        fail = false;
        if(i != generic_offset) {
            for(size_t j=0; j<both_sigs.size(); j++) {
                if(both_sigs[j].is_float) {
                    if(both_sigs[j].size == 8) {
                        if (both_sigs[j].is_complex) {
                            fail = ccompare((double *) test_data[generic_offset][j], (double *) test_data[i][j], vlen, tol_f, absolute_mode);
                        } else {
                            fail = fcompare((double *) test_data[generic_offset][j], (double *) test_data[i][j], vlen, tol_f, absolute_mode);
                        }
                    } else {
                        if (both_sigs[j].is_complex) {
                            fail = ccompare((float *) test_data[generic_offset][j], (float *) test_data[i][j], vlen, tol_f, absolute_mode);
                        } else {
                            fail = fcompare((float *) test_data[generic_offset][j], (float *) test_data[i][j], vlen, tol_f, absolute_mode);
                        }
                    }
                } else {
                    //i could replace this whole switch statement with a memcmp if i wasn't interested in printing the outputs where they differ
                    switch(both_sigs[j].size) {
                        case 8:
                            if(both_sigs[j].is_signed) {
                                fail = icompare((int64_t *) test_data[generic_offset][j], (int64_t *) test_data[i][j], vlen*(both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                            } else {
                                fail = icompare((uint64_t *) test_data[generic_offset][j], (uint64_t *) test_data[i][j], vlen*(both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                            }
                            break;
                        case 4:
                            if(both_sigs[j].is_complex) {
                                if(both_sigs[j].is_signed) {
                                    fail = icompare((int16_t *) test_data[generic_offset][j], (int16_t *) test_data[i][j], vlen*(both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                                } else {
                                    fail = icompare((uint16_t *) test_data[generic_offset][j], (uint16_t *) test_data[i][j], vlen*(both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                                }
                            }
                            else {
                                if (both_sigs[j].is_signed) {
                                    fail = icompare((int32_t *) test_data[generic_offset][j], (int32_t *) test_data[i][j],
                                                    vlen * (both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                                } else {
                                    fail = icompare((uint32_t *) test_data[generic_offset][j], (uint32_t *) test_data[i][j],
                                                    vlen * (both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                                }
                            }
                            break;
                        case 2:
                            if(both_sigs[j].is_signed) {
                                fail = icompare((int16_t *) test_data[generic_offset][j], (int16_t *) test_data[i][j], vlen*(both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                            } else {
                                fail = icompare((uint16_t *) test_data[generic_offset][j], (uint16_t *) test_data[i][j], vlen*(both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                            }
                            break;
                        case 1:
                            if(both_sigs[j].is_signed) {
                                fail = icompare((int8_t *) test_data[generic_offset][j], (int8_t *) test_data[i][j], vlen*(both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                            } else {
                                fail = icompare((uint8_t *) test_data[generic_offset][j], (uint8_t *) test_data[i][j], vlen*(both_sigs[j].is_complex ? 2 : 1), tol_i, absolute_mode);
                            }
                            break;
                        default:
                            fail=1;
                    }
                }
                if(fail) {
                    volk_test_time_t *result = &results->back().results[arch_list[i]];
                    result->pass = !fail;
                    fail_global = true;
                    std::cout << name << ": fail on arch " << arch_list[i] << std::endl;
                }
            }
        }
        arch_results.push_back(!fail);
    }

    double best_time_a = std::numeric_limits<double>::max();
    double best_time_u = std::numeric_limits<double>::max();
    std::string best_arch_a = "generic";
    std::string best_arch_u = "generic";
    for(size_t i=0; i < arch_list.size(); i++)
    {
        if((profile_times[i] < best_time_u) && arch_results[i] && desc.impl_alignment[i] == 0)
        {
            best_time_u = profile_times[i];
            best_arch_u = arch_list[i];
        }
        if((profile_times[i] < best_time_a) && arch_results[i])
        {
            best_time_a = profile_times[i];
            best_arch_a = arch_list[i];
        }
    }

    std::cout << "Best aligned arch: " << best_arch_a << std::endl;
    std::cout << "Best unaligned arch: " << best_arch_u << std::endl;

    if(puppet_master_name == "NULL") {
        results->back().config_name = name;
    } else {
        results->back().config_name = puppet_master_name;
    }
    results->back().best_arch_a = best_arch_a;
    results->back().best_arch_u = best_arch_u;

    return fail_global;
}


// ########################## VOLK PROFILE ############################################

namespace fs = boost::filesystem;

void write_results(const std::vector<volk_test_results_t> *results, const std::string path)
{
    const fs::path config_path(path);
    if (! fs::exists(config_path.parent_path()))
    {
        __android_log_write(ANDROID_LOG_INFO, "volk", "folder doesn't exist, creating");
        std::cout << "Creating " << config_path.parent_path() << "..." << std::endl;
        fs::create_directories(config_path.parent_path());
    }

    std::ofstream config;
    std::cout << "Writing " << path << "..." << std::endl;
    config.open(path.c_str());
    if (!config.is_open()) { //either we don't have write access or we don't have the dir yet
        std::cout << "Error opening file " << path << std::endl;
        __android_log_write(ANDROID_LOG_INFO, "volk", "cannot write to file");
    }

    config << "\
#this file is generated by volk_profile.\n\
#the function name is followed by the preferred architecture.\n\
";

    std::vector<volk_test_results_t>::const_iterator profile_results;
    for(profile_results = results->begin(); profile_results != results->end(); ++profile_results) {
        config << profile_results->config_name << " "
               << profile_results->best_arch_a << " "
               << profile_results->best_arch_u << std::endl;
    }
    config.close();
    __android_log_write(ANDROID_LOG_INFO, "volk", "done writing result");
}

void write_json(std::ofstream &json_file, std::vector<volk_test_results_t> results)
{
    json_file << "{" << std::endl;
    json_file << " \"volk_tests\": [" << std::endl;
    size_t len = results.size();
    size_t i = 0;
    std::vector<volk_test_results_t>::iterator result;
    for(result = results.begin(); result != results.end(); ++result) {
        json_file << "  {" << std::endl;
        json_file << "   \"name\": \"" << result->name << "\"," << std::endl;
        json_file << "   \"vlen\": " << (int)(result->vlen) << "," << std::endl;
        json_file << "   \"iter\": " << result->iter << "," << std::endl;
        json_file << "   \"best_arch_a\": \"" << result->best_arch_a
                  << "\"," << std::endl;
        json_file << "   \"best_arch_u\": \"" << result->best_arch_u
                  << "\"," << std::endl;
        json_file << "   \"results\": {" << std::endl;
        size_t results_len = result->results.size();
        size_t ri = 0;

        std::map<std::string, volk_test_time_t>::iterator kernel_time_pair;
        for(kernel_time_pair = result->results.begin(); kernel_time_pair != result->results.end(); ++kernel_time_pair) {
            volk_test_time_t time = kernel_time_pair->second;
            json_file << "    \"" << time.name << "\": {" << std::endl;
            json_file << "     \"name\": \"" << time.name << "\"," << std::endl;
            json_file << "     \"time\": " << time.time << "," << std::endl;
            json_file << "     \"units\": \"" << time.units << "\"" << std::endl;
            json_file << "    }" ;
            if(ri+1 != results_len) {
                json_file << ",";
            }
            json_file << std::endl;
            ri++;
        }
        json_file << "   }" << std::endl;
        json_file << "  }";
        if(i+1 != len) {
            json_file << ",";
        }
        json_file << std::endl;
        i++;
    }
    json_file << " ]" << std::endl;
    json_file << "}" << std::endl;
}

volk_test_params_t test_params(1e-6f, 327.f, 131071, 1987, true, "");

extern "C"
JNIEXPORT jstring JNICALL
Java_net_bastibl_volk_MainActivity_runVolk(JNIEnv *env, jobject thiz) {

    std::string config_file = getenv("EXTERNAL_STORAGE");
    std::ofstream json_file;
    json_file.open(config_file + "/volk/volk_results.json");
    config_file += "/volk/volk_config";


    // Run tests
    std::vector<volk_test_results_t> results;

    // Initialize the list of tests
    std::vector<volk_test_case_t> test_cases = init_test_list(test_params);

    for(unsigned int ii = 0; ii < test_cases.size(); ++ii) {

        volk_test_case_t test_case = test_cases[ii];

        try {
            __android_log_write(ANDROID_LOG_INFO, "volk", test_case.name().c_str());
            run_volk_tests(test_case.desc(), test_case.kernel_ptr(), test_case.name(),
                           test_case.test_parameters(), &results, test_case.puppet_master_name());
        }
        catch (std::string &error) {
            std::cerr << "Caught Exception in 'run_volk_tests': " << error << std::endl;
        }
    }

    write_json(json_file, results);
    json_file.close();

    write_results(&results, config_file);

    return env->NewStringUTF("done");
}