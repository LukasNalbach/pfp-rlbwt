/**
 * part of LukasNalbach/pfp-rlbwt
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// explicit instantiations of pfp_rlbwt::build and of the text sources for one symbol type; this file is compiled
// once per symbol width (PFP_RLBWT_SYMBOL_BYTES in [1, 8], see CMakeLists.txt)

#include <pfp_rlbwt/pfp_rlbwt.hpp>

namespace pfp_rlbwt {

using instantiated_symbol_t = uint_of_width_t<PFP_RLBWT_SYMBOL_BYTES>;

template class file_source<instantiated_symbol_t>;
template class memory_source<instantiated_symbol_t>;
template result build<instantiated_symbol_t>(const text_source<instantiated_symbol_t>&, const std::string&,
    const params&);

} // namespace pfp_rlbwt
