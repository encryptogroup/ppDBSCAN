//
// Created by helenmoellering on 02.03.20. Thanks to Lennart Braun and Amos Treiber for the code.
//

#ifndef ABY_CIRCUITWRAPPER_H
#define ABY_CIRCUITWRAPPER_H

#include <algorithm>
#include <bitset>
#include <cassert>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <type_traits>
#include <abycore/circuit/booleancircuits.h>
#include <abycore/circuit/share.h>

using bytes_t = std::vector<uint8_t>;
using share_p = std::shared_ptr<share>;

class CircuitWrapper
{
public:
  CircuitWrapper(Circuit *circ) : circ_(circ) {}

  template <size_t len>
  std::vector<share_p> PutINGates(std::bitset<len> bits, e_role role)
  {
    std::vector<share_p> shares(len);
    for (size_t i{0}; i < len; ++i)
    {
      shares[i] = share_p(circ_->PutINGate(static_cast<uint8_t>(bits[i]), 1, role));
    }
    return shares;
  }

  std::vector<share_p> PutINGates(const bytes_t &in, e_role role)
  {
    return PutINGates(in, 8, role);
  }

  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  std::vector<share_p> PutINGates(const std::vector<T> &input, size_t size, e_role role)
  {
    std::vector<share_p> output(input.size());
    std::transform(input.begin(), input.end(), output.begin(),
                   [this, size, role] (uint8_t in) -> share_p
                   { return this->PutINGate(in, static_cast<uint32_t>(size), role); });
    return output;
  }

  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  share_p PutINGate(T input, size_t size, e_role role)
  {
    assert(size <= sizeof(input) * 8);
    assert(size <= std::numeric_limits<uint32_t>::max());
    return share_p(circ_->PutINGate(input, static_cast<uint32_t>(size), role));
  }

  share_p PutDummyINGate(size_t size)
  {
    assert(size <= std::numeric_limits<uint32_t>::max());
    return share_p(circ_->PutDummyINGate(static_cast<uint32_t>(size)));
  }

  std::vector<share_p> PutDummyINGates(size_t num, size_t size)
  {
    assert(size <= std::numeric_limits<uint32_t>::max());
    std::vector<share_p> shares(num);
    for (size_t i{0}; i < num; ++i)
    {
      shares[i] = share_p(circ_->PutDummyINGate(static_cast<uint32_t>(size)));
    }
    return shares;
  }

  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  share_p PutSharedSIMDINGate(uint32_t num, T* invals, uint32_t bitlen) {
    return share_p (circ_->PutSharedSIMDINGate(num, invals, bitlen));
  }

  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  share_p PutSIMDINGate(uint32_t num, T* invals, size_t bitlen, e_role role) {
    return share_p (circ_->PutSIMDINGate(num, invals, bitlen,role));
  }


  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  share_p PutSIMDCONSGate(uint32_t val, T input, size_t size)
  {
    return share_p(circ_->PutSIMDCONSGate(val,input, static_cast<uint32_t>(size)));
  }

  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  share_p PutSharedINGate(T val, uint32_t bitlen)
  {
    return share_p(circ_->PutSharedINGate(val, bitlen));
  }

  share_p PutOUTGate(share_p ina, e_role role)
  {
    return share_p(circ_->PutOUTGate(ina.get(), role));
  }

  share_p PutSharedOUTGate(share_p ina)
  {
    return share_p(circ_->PutSharedOUTGate(ina.get()));
  }

  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  share_p PutCONSGate(T input, size_t size)
  {
    //assert(size <= sizeof(input) * 8);
    //assert(size <= std::numeric_limits<uint32_t>::max());
    return share_p(circ_->PutCONSGate(input, static_cast<uint32_t>(size)));
  }

  share_p PutANDGate(share_p ina, share_p inb)
  {
    return share_p(circ_->PutANDGate(ina.get(), inb.get()));
  }

  share_p PutXORGate(share_p ina, share_p inb)
  {
    return share_p(circ_->PutXORGate(ina.get(), inb.get()));
  }

  share_p PutORGate(share_p ina, share_p inb)
  {
    return share_p(((BooleanCircuit*)circ_)->PutORGate(ina.get(), inb.get()));
  }

  share_p PutADDGate(share_p ina, share_p inb)
  {
    return share_p(circ_->PutADDGate(ina.get(), inb.get()));
  }

  share_p PutSUBGate(share_p ina, share_p inb)
  {
    return share_p(circ_->PutSUBGate(ina.get(), inb.get()));
  }

  share_p PutMULGate(share_p ina, share_p inb)
  {
    return share_p(circ_->PutMULGate(ina.get(), inb.get()));
  }

  share_p PutGTGate(share_p ina, share_p inb)
  {
    return share_p(circ_->PutGTGate(ina.get(), inb.get()));
  }

  share_p PutEQGate(share_p ina, share_p inb)
  {
    return share_p(circ_->PutEQGate(ina.get(), inb.get()));
  }

  share_p PutMUXGate(share_p ina, share_p inb, share_p sel)
  {
    return share_p(circ_->PutMUXGate(ina.get(), inb.get(), sel.get()));
  }

  share_p PutHammingWeightGateRec(uint32_t *wires,size_t size){
    return share_p(((BooleanCircuit*)circ_)->PutHammingWeightGateRec(wires,size));
  }

  share_p PutHammingWeightGate(share_p in) {
    return share_p(((BooleanCircuit*)circ_)->PutHammingWeightGate(in.get()));
  }

  share_p PutSubsetGate(share_p in, uint32_t* pos, uint32_t nvals_out) {
    return share_p(circ_->PutSubsetGate(in.get(), pos, nvals_out));
  }

  share_p PutSplitterGate(share_p in) {
    return share_p(circ_->PutSplitterGate(in.get()));
  }

  share_p PutRepeaterGate(uint32_t val, share_p in) {
    return share_p(circ_->PutRepeaterGate(val,in.get()));
  }

  share_p PutCombinerGate(share_p in) {
    return share_p(circ_->PutCombinerGate(in.get()));
  }

  uint32_t PutCombinerGate(std::vector<uint32_t> in) {
    return circ_->PutCombinerGate(in);
  }

  //share_p PutMinGate(share_p in, uint32_t nvals_out) {
    //return share_p(circ_->PutMinGate(in.get(), nvals_out));
  //}

  share_p PutROTLGate(share_p in, size_t k)
  {
    size_t bitlength{in->get_bitlength()};
    auto wires{in->get_wires()};
    k %= bitlength;
    auto new_first{std::prev(wires.end(), static_cast<ptrdiff_t>(k))};
    std::rotate(wires.begin(), new_first, wires.end());
    return std::make_shared<boolshare>(wires, circ_);
  }

  share_p PutReduceANDGates(std::vector<share_p> ins)
  {
    return PutReduceGates([this] (share_p ina, share_p inb) -> share_p
                          { return this->PutANDGate(ina, inb); }, ins);
  }

  share_p PutReduceGates(std::function<share_p(share_p, share_p)> PutGate,
                         std::vector<share_p> ins)
  {
    assert(!ins.empty());
    auto result{std::accumulate(ins.cbegin() + 1, ins.cend(), ins[0],
                                [PutGate] (share_p acc, share_p val) -> share_p
                                { return PutGate(acc, val); })};
    return result;
  }


  // debug gates
  share_p PutPrintValueGate(share_p in, std::string infostring)
  {
    return share_p(circ_->PutPrintValueGate(in.get(), infostring));
  }

  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  share_p PutAssertGate(share_p in, T value, size_t size)
  {
    return share_p(circ_->PutAssertGate(in, value, size));
  }

  share_p PutA2YGate(share_p in) {
    return share_p(circ_->PutA2YGate(in.get()));
  }

  share_p PutB2YGate(share_p in) {
    return share_p(circ_->PutB2YGate(in.get()));
  }

  share_p PutA2BGate(share_p in, Circuit* yaosharingcircuit) {
    return share_p(circ_->PutA2BGate(in.get(),yaosharingcircuit));
  }

  share_p PutB2AGate(share_p in) {
    return share_p(circ_->PutB2AGate(in.get()));
  }

  share_p PutY2AGate(share_p in, Circuit* boolsharingcircuit) {
    return share_p(circ_->PutY2AGate(in.get(),boolsharingcircuit));
  }

  share_p PutY2BGate(share_p in) {
    return share_p(circ_->PutY2BGate(in.get()));
  }

  static share_p get_wire_ids_as_share(share_p share, size_t index)
  {
    assert(index <= std::numeric_limits<uint32_t>::max());
    return share_p(share->get_wire_ids_as_share(static_cast<uint32_t>(index)));
  }

  // share conversion
  std::vector<share_p> split_share(share_p in, size_t bits_per_share)
  {
    assert(in->get_bitlength() % bits_per_share == 0);
    size_t num_shares{in->get_bitlength() / bits_per_share};
    std::vector<share_p> out(num_shares);
    auto in_wires{in->get_wires()};
    for (size_t i{0}; i < num_shares; ++i)
    {
      std::vector<uint32_t> wires(in_wires.cbegin() + static_cast<ptrdiff_t>(i*bits_per_share),
                                  in_wires.cbegin() + static_cast<ptrdiff_t>((i+1)*bits_per_share));

      out[i] = std::make_shared<boolshare>(wires, circ_);
    }
    return out;
  }

  share_p combine_shares(std::vector<share_p> in)
  {
    size_t num_bits{std::accumulate(in.cbegin(), in.cend(), size_t{0},
                                    [] (size_t acc, share_p s) -> size_t
                                    { return acc + s->get_bitlength(); })};
    assert(num_bits <= 64);
    std::vector<uint32_t> wires;
    for (auto &share : in)
    {
      auto in_wires{share->get_wires()};
      wires.insert(wires.cend(), in_wires.cbegin(), in_wires.cend());
    }
    return std::make_shared<boolshare>(wires, circ_);
  }

  template <typename T,
      typename = std::enable_if_t<std::is_unsigned<T>::value>>
  share_p replaceValueInSIMD(share_p shr,uint32_t pos,T new_value)
  {

    //return share_p(circ_->replaceValueInSIMD(shr.get(),pos,new_value.get()));
    /*std::vector<uint32_t> wireIds = shr.get()->get_wires();

    for (uint32_t i = 0; i < circ_.m_vGates[wireIds[0]].sharebitlen; i++) {
      GATE *gate = &(m_vGates[gateid]);
      //memcpy(((uint8_t *)gate->gs.val) + pos * sharebytelen, new_value,
      //       inbytelen);
    }*/
  }


  std::vector<share_p> bytes_to_shares(const bytes_t &in)
  {
    std::vector<share_p> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [this] (uint8_t b) -> share_p
                   { return this->PutCONSGate(b, 8); });
    return out;
  }

  bytes_t shares_to_bytes(const std::vector<share_p> &in)
  {
    bytes_t out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [] (share_p b) -> uint8_t
                   {
                     assert(b->get_bitlength() == 8);
                     return b->template get_clear_value<uint8_t>();
                   });
    return out;
  }
// private:
  Circuit* circ_;
};

using CircuitW_p = std::shared_ptr<CircuitWrapper>;

#endif // ABY_CIRCUITWRAPPER_H
