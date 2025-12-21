//
// Created by Noah bittermann on 12/7/25.
//

#ifndef STENCIL_H
#define STENCIL_H

#include <numeric>
#include <cassert>

#include "pragma_shorts.h"

class SiteNeighbors {
public:
    SiteNeighbors() = default;
    explicit SiteNeighbors(int d) {
        p.resize(d);
        m.resize(d);
    }

    std::vector<int> p;
    std::vector<int> m;
};

template<int d>
class Stencil {
public:
    std::vector<SiteNeighbors> neighbors;

    explicit Stencil(const std::array<int,d>& dims) :
    _dims(dims),
    _strides(d) {

        int sites = std::accumulate(begin(dims), end(dims), 1, std::multiplies<>());

        _strides[d-1] = 1; for (int a = d - 2; a >= 0; --a) _strides[a] = _strides[a+1]*_dims[a+1];

        neighbors.reserve(sites);
        for(int ss=0; ss < sites; ++ss) neighbors.emplace_back(d);

        thread_for(ss, sites, {
            std::vector<int> coords = ravel(ss);
            for(int dir=0; dir<d; ++dir) {
                std::vector<int> plus_neigh = coords;
                std::vector<int> minus_neigh = coords;

                plus_neigh[dir] = (plus_neigh[dir] + 1 == _dims[dir]) ? 0 : plus_neigh[dir] + 1;
                minus_neigh[dir] = (minus_neigh[dir] - 1 == -1) ? _dims[dir] - 1 : minus_neigh[dir] - 1;

                neighbors[ss].p[dir] = unravel(plus_neigh);
                neighbors[ss].m[dir] = unravel(minus_neigh);
            }
        });
    }

    int unravel(const std::vector<int>& coords) const {
        assert(coords.size() == d);
        int ss=0; for (int a=0; a<d; ++a) ss += coords[a]*_strides[a];
        return ss;
    }

    std::vector<int> ravel(int ss) const {
        std::vector<int> coords;
        int rem = ss;
        for (int a = 0; a < d; ++a) {
            coords.push_back(rem / _strides[a]);
            rem %= _strides[a];
        }
        return coords;
    }

private:
    std::array<int, d> _dims;
    std::vector<int> _strides;

};

#endif //STENCIL_H

