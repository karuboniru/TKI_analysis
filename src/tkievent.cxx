#include "tkievent.h"
#include <stdexcept>

size_t TKIEvent::count_in(int id) const { return in.count(id); }
size_t TKIEvent::count_out(int id) const { return out.count(id); }
size_t TKIEvent::count_post(int id) const { return post.count(id); }

void TKIEvent::add_in(int id, const TLorentzVector particle) {
  in.insert({id, particle});
}
void TKIEvent::add_out(int id, const TLorentzVector particle) {
  out.insert({id, particle});
}
void TKIEvent::add_post(int id, const TLorentzVector particle) {
  post.insert({id, particle});
  ids_post.insert(id);
}
const std::set<int> &TKIEvent::get_ids_post() const { return ids_post; }

const TLorentzVector &TKIEvent::get_leading(int id) const {
  const TLorentzVector *leading = nullptr;
  for (const auto &p : post_range(id)) {
    if (leading == nullptr || p.second.P() > leading->P()) {
      leading = &p.second;
    }
  }
  if (leading == nullptr) {
    throw std::runtime_error("No leading particle found");
  }
  return *leading;
}