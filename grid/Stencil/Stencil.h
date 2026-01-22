#ifndef STENCIL_H
#define STENCIL_H

#include <array>
#include <vector>
#include <numeric>
#include <functional>
#include <cassert>

#include "../openmp/pragma_shorts.h"
#include "../simd/simd_types.h"

template<int d>
struct BlockNeighbors {
    std::array<int, d> p;
    std::array<int, d> m;
};

template<int d>
class BlockStencil {
public:
    std::vector<BlockNeighbors<d>> neighbors;

    explicit BlockStencil(const std::array<int, d>& dims_sites) :
    _dims_sites(dims_sites),
    _dims_blocks(dims_sites),
    _strides_blocks(d)
    {
        constexpr int pack_dim = 0;
        constexpr int W = simd::W;

        if ((_dims_sites[pack_dim] % W) != 0) {
            throw std::runtime_error("BlockStencil: dims[0] must be a multiple of simd::W when packing across x.");
        }

        _dims_blocks[pack_dim] = _dims_sites[pack_dim] / W;

        const int blocks = std::accumulate(
                    _dims_blocks.begin(), _dims_blocks.end(),
                    1, std::multiplies<>()
        );

        _strides_blocks[0] = 1;
        for (int a = 1; a < d; ++a) {
            _strides_blocks[a] = _strides_blocks[a - 1] * _dims_blocks[a - 1];
        }

        neighbors.resize(blocks);

        thread_for(b, blocks, {
            auto coords = ravel_block(b);
            for (int dir = 0; dir < d; ++dir) {
                auto plus  = coords;
                auto minus = coords;

                plus[dir]  = (plus[dir]  + 1 == _dims_blocks[dir]) ? 0 : plus[dir]  + 1;
                minus[dir] = (minus[dir] - 1 <  0)                 ? _dims_blocks[dir] - 1 : minus[dir] - 1;

                neighbors[b].p[dir] = unravel_block(plus);
                neighbors[b].m[dir] = unravel_block(minus);
            }
        });
    }

    int unravel_block(const std::array<int,d>& coords) const {
        int b = 0;
        for (int a = 0; a < d; ++a) b += coords[a] * _strides_blocks[a];
        return b;
    }

    std::array<int,d> ravel_block(int b) const {
        std::array<int,d> coords{};
        int rem = b;
        for (int a = d - 1; a >= 0; --a) {
            coords[a] = rem / _strides_blocks[a];
            rem %= _strides_blocks[a];
        }
        return coords;
    }

    const std::array<int,d>& dims_sites()  const { return _dims_sites;  }
    const std::array<int,d>& dims_blocks() const { return _dims_blocks; }


private:
    std::array<int,d> _dims_sites;
    std::array<int,d> _dims_blocks;
    std::vector<int>  _strides_blocks;
};

#endif //STENCIL_H

