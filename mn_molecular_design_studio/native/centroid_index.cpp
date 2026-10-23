// Exact, ordered sphere exclusion on packed fingerprints.
// OpenVS's best-first rule is retained: the earliest selected representative
// within the threshold owns each compound. Prefix filtering is lossless.
#include <algorithm>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

namespace {
using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u128 = unsigned __int128;
constexpr u64 missing = std::numeric_limits<u64>::max();
thread_local std::string error_message;
using Progress = int (*)(u64, u64, u64);

struct Scratch {
  std::vector<u32> seen;
  std::vector<u32> bits;
  u32 epoch = 0;
};

void prefix(const u64 *fp, u32 words, u32 count, u64 p, u64 q,
            const std::vector<u32> &rank, std::vector<u32> &bits) {
  bits.clear();
  for (u32 w = 0; w < words; ++w) {
    u64 value = fp[w];
    while (value) {
      bits.push_back(w * 64 + __builtin_ctzll(value));
      value &= value - 1;
    }
  }
  std::sort(bits.begin(), bits.end(), [&rank](u32 a, u32 b) { return rank[a] < rank[b]; });
  if (count) {
    const u32 required = static_cast<u32>((u128(p) * count + q - 1) / q);
    bits.resize(count - required + 1);
  }
}

bool similar(const u64 *a, const u64 *b, u32 words, u32 ac, u32 bc, u64 p, u64 q, bool zero_matches) {
  if (!ac && !bc) return zero_matches;
  u32 intersection = 0;
  for (u32 w = 0; w < words; ++w) intersection += __builtin_popcountll(a[w] & b[w]);
  return u128(intersection) * (u128(p) + q) >= u128(p) * (u64(ac) + bc);
}

template <typename Index>
int cluster(const u64 *fps, u64 n, u32 words, const u64 *order,
            const u64 *weights,
            u64 p, u64 q, u32 threads, u32 block_size, u32 zero_matches,
            u64 *representatives, u64 *sizes, u64 *summary, Progress progress) {
  const u32 nbits = words * 64;
  std::vector<u32> popcounts(n);
  std::vector<std::vector<u64>> frequencies(threads, std::vector<u64>(nbits, 0));
  #pragma omp parallel num_threads(threads)
  {
    auto &freq = frequencies[omp_get_thread_num()];
    #pragma omp for schedule(static)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i) {
      u32 count = 0;
      for (u32 w = 0; w < words; ++w) {
        u64 value = fps[u64(i) * words + w];
        count += __builtin_popcountll(value);
        while (value) {
          ++freq[w * 64 + __builtin_ctzll(value)];
          value &= value - 1;
        }
      }
      popcounts[i] = count;
    }
  }
  std::vector<u64> frequency(nbits, 0);
  for (const auto &freq : frequencies)
    for (u32 bit = 0; bit < nbits; ++bit) frequency[bit] += freq[bit];
  frequencies.clear();
  std::vector<u32> bit_order(nbits), rank(nbits);
  std::iota(bit_order.begin(), bit_order.end(), 0);
  std::sort(bit_order.begin(), bit_order.end(), [&frequency](u32 a, u32 b) {
    return frequency[a] != frequency[b] ? frequency[a] < frequency[b] : a < b;
  });
  for (u32 i = 0; i < nbits; ++i) rank[bit_order[i]] = i;
  std::vector<std::vector<Index>> postings(nbits);
  std::vector<Scratch> scratch(threads);
  for (auto &state : scratch) state.bits.reserve(nbits);
  std::vector<u32> new_bits;
  new_bits.reserve(nbits);
  std::vector<u64> matches(block_size, missing);
  u64 k = 0, checks = 0, posting_count = 0, zero_rep = missing;
  if (progress && progress(0, 0, 0)) return 1;

  for (u64 start = 0; start < n; start += block_size) {
    const u64 end = std::min(n, start + block_size), old_k = k;
    for (auto &state : scratch) state.seen.resize(old_k, 0);
    u64 block_checks = 0;
    // The index is immutable for this phase. Each query can run independently.
    #pragma omp parallel for schedule(dynamic, 16) num_threads(threads) reduction(+:block_checks)
    for (std::int64_t offset = 0; offset < static_cast<std::int64_t>(end - start); ++offset) {
      const u64 source = order[start + offset];
      const auto *fp = fps + source * words;
      const u32 count = popcounts[source];
      u64 best = missing;
      if (!count) {
        best = zero_matches ? zero_rep : missing;
      } else if (old_k) {
        auto &state = scratch[omp_get_thread_num()];
        if (++state.epoch == 0) {
          std::fill(state.seen.begin(), state.seen.end(), 0);
          ++state.epoch;
        }
        prefix(fp, words, count, p, q, rank, state.bits);
        const u64 minimum = static_cast<u64>((u128(p) * count + q - 1) / q);
        const u64 maximum = static_cast<u64>(std::min(u128(nbits), u128(q) * count / p));
        for (const u32 bit : state.bits) {
          for (const Index raw : postings[bit]) {
            const u64 candidate = raw;
            if (candidate >= best) break;
            if (state.seen[candidate] == state.epoch) continue;
            state.seen[candidate] = state.epoch;
            const u64 other = representatives[candidate];
            const u32 other_count = popcounts[other];
            if (other_count < minimum || other_count > maximum) continue;
            ++block_checks;
            if (similar(fp, fps + other * words, words, count, other_count, p, q, zero_matches)) best = candidate;
          }
        }
      }
      matches[offset] = best;
    }
    checks += block_checks;
    // Resolve newly created representatives in input order. This preserves
    // exactly the result of a serial best-first pass; chunks are never merged.
    for (u64 position = start; position < end; ++position) {
      const u64 source = order[position];
      const u32 count = popcounts[source];
      u64 match = matches[position - start];
      if (match == missing) {
        for (u64 candidate = old_k; candidate < k; ++candidate) {
          const u64 other = representatives[candidate];
          const u32 other_count = popcounts[other];
          if (u128(p) * count > u128(q) * other_count || u128(p) * other_count > u128(q) * count) continue;
          ++checks;
          if (similar(fps + source * words, fps + other * words, words, count, other_count, p, q, zero_matches)) {
            match = candidate;
            break;
          }
        }
      }
      if (match == missing) {
        match = k++;
        representatives[match] = source;
        sizes[match] = 0;
        if (!count && zero_rep == missing) zero_rep = match;
      }
      sizes[match] += weights ? weights[source] : 1;
    }
    for (u64 candidate = old_k; candidate < k; ++candidate) {
      const u64 source = representatives[candidate];
      prefix(fps + source * words, words, popcounts[source], p, q, rank, new_bits);
      for (const u32 bit : new_bits) postings[bit].push_back(static_cast<Index>(candidate));
      posting_count += new_bits.size();
    }
    summary[0] = k;
    summary[1] = checks;
    summary[2] = posting_count;
    summary[3] = threads;
    if (progress && progress(end, k, checks)) return 1;
  }
  return 0;
}
}  // namespace

extern "C" const char *studio_centroid_error() { return error_message.c_str(); }

extern "C" int studio_centroid_cluster(const u64 *fps, u64 n, u32 words, const u64 *order,
                                      const u64 *weights,
                                      u64 p, u64 q, u32 threads, u32 block_size, u32 zero_matches,
                                      u64 *representatives, u64 *sizes, u64 *summary, Progress progress) {
  try {
    if (!n || !words || !threads || !block_size || !p || p > q)
      throw std::runtime_error("Invalid centroid backend arguments.");
    for (u64 i = 0; i < n; ++i)
      if (order[i] >= n) throw std::runtime_error("Centroid ordering references an invalid row.");
    if (n <= std::numeric_limits<u32>::max())
      return cluster<u32>(fps, n, words, order, weights, p, q, threads, block_size, zero_matches, representatives, sizes, summary, progress);
    return cluster<u64>(fps, n, words, order, weights, p, q, threads, block_size, zero_matches, representatives, sizes, summary, progress);
  } catch (const std::exception &exc) {
    error_message = exc.what();
    return -1;
  }
}
