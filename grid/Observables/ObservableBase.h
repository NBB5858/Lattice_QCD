#ifndef OBSERVABLEBASE_H
#define OBSERVABLEBASE_H

class ObservableBase {
public:

    explicit ObservableBase(const std::size_t nobs) : _nobs(nobs) {_cache.reserve(_nobs);}
    virtual ~ObservableBase() = default;

    const std::vector<float>& cache() const { return _cache; }

    float mean(float burnin) const {
        if (_cache.empty()) return 0.0;

        const std::size_t n = _cache.size();
        const auto start = static_cast<std::size_t>(std::ceil(burnin * n));

        float s = 0.0f;
        for (std::size_t i = start; i < n; ++i) s += _cache[i];

        const std::size_t denom = (n > start) ? (n - start) : 1;
        return s / static_cast<float>(denom);
    }

    float abs_mean(float burnin) const {
        if (_cache.empty()) return 0.0;

        const std::size_t n = _cache.size();
        const auto start = static_cast<std::size_t>(std::ceil(burnin * n));

        float s = 0.0f;
        for (std::size_t i = start; i < n; ++i) {
            s += std::abs(_cache[i]);
        }

        const std::size_t denom = (n > start) ? (n - start) : 1;
        return s / static_cast<float>(denom);
    }

    std::vector<float> _cache;
protected:
    std::size_t _nobs;
};

#endif //OBSERVABLEBASE_H