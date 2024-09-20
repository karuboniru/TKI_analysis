#pragma once

#include <TLorentzVector.h>
#include <set>
#include <unordered_map>
#include <utility>

template <typename U> class equal_range_iterable {
public:
  template <typename T>
  equal_range_iterable(T &&m, int id) : range(m.equal_range(id)) {}

  auto &&begin() { return range.first; }

  auto &&end() { return range.second; }

private:
  decltype(std::declval<U>().equal_range(std::declval<int>())) range;
  int id;
};

template <typename U>
equal_range_iterable(U &&, int) -> equal_range_iterable<U>;

class NeutrinoEvent {
public:
  size_t count_in(int) const;
  size_t count_out(int) const;
  size_t count_post(int) const;

  void add_in(int, const TLorentzVector);
  void add_out(int, const TLorentzVector);
  void add_post(int, const TLorentzVector);

  auto in_range(int id) const { return equal_range_iterable(in, id); }
  auto out_range(int id) const { return equal_range_iterable(out, id); }
  auto post_range(int id) const { return equal_range_iterable(post, id); }

  const std::set<int> &get_ids_post() const;

  const TLorentzVector & get_leading(int id) const;

private:
  std::unordered_multimap<int, TLorentzVector> in, out, post;
  std::set<int> ids_post;
};