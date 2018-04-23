/* clang-format off */
/*
iji2dgrid : IncredibleJunior Implicit 2D Grid

iji2dgrid is a grid that can be used to accelerate broadphase queries, or for
coarse collision detection.

The grid is optimized for (low) memory consumption first and foremost where our
use case consisted of (relative) low number of columns/rows (<20) and large
number of items (>16000).

The memory requirements can be calculated as follows (64-bit version):

bytes_needed =
   (num_rows + num_columns) * ( ((max_num_items+63) / 64) * 8) (for bitarrays)
   + max_num_items*max(4, item_size) (for per item data)
   + 7 (bytes to ensure natural alignment of bitarrays)

The grid has a user configurable per grid-item size and optional per cell data
and only uses memory that it is provided during initialization. Whose size depends
on the result of the configuration of the grid, that is the number of items,
size in bytes of each item, per cell data (if any) and how many rows and columns
that make up the grid.

When an item is added to the grid, a (grid-)handle is constructed that is used to
reference that item for later retrieval of userdata and removal from the grid.

A grid-handle is encoded with item index and which cells the item overlaps
(cell_start_x, cell_start_y, cell_end_x, cell_end_y). How the bit-configuration
of a handle looks like is decided by the grid-configuration. The handle's
bit-configuration is derived by first reserving the number of bits needed to
represent the max number of items then the remaining bits is divided into four
equally large number of bits, where we store the cell coordinates the handles item
overlaps. We do this in order to be able to handle the worst case, where a item
overlap all cells (one could have a separate bit configuration for columns and
rows, if the grid is not square, yielding a higher max num items count).

A grid-handle is by default 64-bit (unless IJI2DGRID_32BIT_HANDLES is defined)
but if a grid configuration yields a handle configuration that can be represented
in 32-bits a flag (IJI2DGRID_FF_32BIT_CLEAN_HANDLES) is set and the handle can
safely be cast down to 32-bits without loosing information.

The grid is divided, by uniformly sized cells, into rows and columns.

Given an AABB extent and cell size we divide the grid like:

                                       cell(num_columns-1, num_rows-1)
   cell(0, num_rows-1)          AABB-max
                   .--------------.  <.
                   |              |   |
                   |__ __ ...     |   | height
                   |  |  |        |   |
                   '--------------'  <'
            cell(0,0)              cell(num_columns-1, 0)
            AABB-min
                   ^--------------^ width

Given this we can calculate unique cell id's for the cells by:
   cellid = column + row * num_columns

and reverse by:
   column = cellid % num_columns
   row = cellid / num_columns

Each row and column is represented by a bitarray of the configured max number
of items num bits (rounded up to nearest word size). Adding an item to the grid
will mark that items index in each of its overlapping columns and rows. Making
an overlap query will first mask together (bitwise OR) two bitfields of the
overlapping columns and rows respectively and then by using bitwise AND between
these column and row bitfields gives us the items present in the overlap.
This technique is written about in 'Real-time Collision Detection'
by Christer Ericson (chapter 7.1.5 Implicit Grids) [1]

This file provides both the interface and the implementation.

The grid is implemented as a stb-style header-file library[2]
which means that in *ONE* source file, put:

#define IJI2DGRID_IMPLEMENTATION
// if no external dependencies is wanted (on assert.h, string.h, math.h)
// you can override the following before including (with implementation defined)
#ifndef IJI2DGRID_ceilf      custom_ceilf
#ifndef IJI2DGRID_memset     custom_memset
#define IJI2DGRID_assert     custom_assert
#include "iji2dgrid.h"

Other source files should just include iji2dgrid.h

NB: the implementation uses intrinsics/builtins for Find First Set(ffs) and
Count Leading Zeroes(clz) and is currently only implemented for Windows.

EXAMPLES/UNIT TESTS
Usage examples+tests is at the bottom of the file in the IJI2DGRID_TEST section.

LICENSE
See end of file for license information

REFERENCES

[1] Real-Time Collision Detection, Christer Ericson, ISBN:9780080474144
[2] https://github.com/nothings/stb
*/

#ifndef IJI2DGRID_INCLUDED_H
#define IJI2DGRID_INCLUDED_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(IJI2DGRID_STATIC)
   #define IJI2DGRID_API static
#else
   #define IJI2DGRID_API extern
#endif

#if defined(IJI2DGRID_32BIT_HANDLES)
   typedef unsigned iji2dgrid_handle;
#else
   #ifdef _MSC_VER
      typedef unsigned __int64 iji2dgrid_handle;
   #else
      typedef unsigned long long iji2dgrid_handle;
   #endif
#endif

/* grid-items handles can safely casted to 32-bit without losing information */
#define IJI2DGRID_FF_32BIT_CLEAN_HANDLES (1<<0)

struct iji2dgrid {
   unsigned short ncols;
   unsigned short nrows;
   unsigned char nbits_per_coord;
   unsigned char nbits_items;

   unsigned short feature_flags; /* IJI2DGRID_FF_32BIT_CLEAN_HANDLES */

   float cell_size;
   float gridmin[2];
   float gridmax[2];

   unsigned items_capacity; /* power of two not required/enforced */
   unsigned items_mask;

   unsigned short item_size_bytes;
   unsigned short celldata_size_bytes;

   unsigned ngridwords_per_row_col;
   unsigned data_freelist;

   void *data; /* memory provided during initialization */
   void *celldata;
};

/* NB: num _usable_ items is max_num_items-1 (0 is reserved for invalid handle) */
IJI2DGRID_API int iji2dgrid_valid_configuration(float width, float height, float cell_size,
   unsigned max_num_items, unsigned item_size_bytes, unsigned cellitem_size_bytes, unsigned cellitem_alignment);

/* returns the maximum number of items that can be serviced (where one is reserved for invalid handle, i.e. 0) */
IJI2DGRID_API unsigned iji2dgrid_max_num_items(float width, float height, float cell_size);

/* computes the number of columns and rows needed given extents and cell size */
IJI2DGRID_API void iji2dgrid_dimensions(float width, float height,
   float cell_size, unsigned *ncols, unsigned *nrows);

/* NB: sizeof struct iji2dgrid is _not_ included */
IJI2DGRID_API unsigned iji2dgrid_memory_size_needed(float width, float height, float cell_size,
   unsigned max_num_items, unsigned item_size_bytes,
   unsigned cellitem_size_bytes, unsigned cellitem_alignment);

/*
initialize the grid with a configuration that is assumed to pass
a 'iji2dgrid_valid_configuration' check.

if no per-cell data is needed then pass 0 to cellitem_size_bytes and cellitem_alignment

the grid-items will be located at the start of the provided memory and needs to
be aligned to the required alignment of the their type.
*/
IJI2DGRID_API void iji2dgrid_init(struct iji2dgrid *self,
   float gridminx, float gridminy, float gridmaxx, float gridmaxy,
   float cell_size, unsigned max_num_items, unsigned item_size_bytes,
   unsigned cellitem_size_bytes, unsigned cellitem_alignment,
   void *memory);

/* reset to initial state */
IJI2DGRID_API void iji2dgrid_reset(struct iji2dgrid *self);

/*
returns a valid pointer to a memory area of 'item_size_bytes', which was configured/set
in 'iji2dgrid_init', bytes and a valid handle on success.

returns a pointer to a non-valid valid memory area (i.e. 0/null) and fills the
handle with invalid handle (0) on failure
*/
IJI2DGRID_API void *iji2dgrid_add_box(struct iji2dgrid *self, float minx, float miny,
   float maxx, float maxy, iji2dgrid_handle *handle_out);
/* add circle by its aabb-extents */
#define iji2dgrid_add_cirle(self, centerx, centery, radius, handleout) iji2dgrid_add_box((self), (centerx)-(radius), (centery)-(radius), (centerx)+(radius), (centery)+(radius), (handleout))

/* removes the item from the grid. NB: all non-zero handles is interpreted as a valid handles */
IJI2DGRID_API void iji2dgrid_remove(struct iji2dgrid *self, iji2dgrid_handle handle);

/* returns a grid coordinate (x(i.e column), y(i.e. row)) */
IJI2DGRID_API unsigned iji2dgrid_coordinate(struct iji2dgrid *self, float posx, float posy);

/*
returns the number of overlaps and, if 'max_num_handles' is non-zero, fills
handles with pointers to grid-items, which is user-defined and has to be
cast to the correct type before use.

if storage for 'handles' is to be preallocated or only a overlap count is
needed pass 0 for 'handles' and 'max_num_handles'.

NB: the return value may be larger than the maximum number of handles
written(which is bound by 'max_num_handles') so be sure to clamp the range
when iterating over the result.
*/
IJI2DGRID_API unsigned iji2dgrid_get_overlaps(struct iji2dgrid *self,
   unsigned coord_start, unsigned coord_end,
   void **handles, unsigned max_num_handles);

/* helpers for overlaps with boxes and circles (by its aabb-extents) */
#define iji2dgrid_get_box_overlaps(self, minx, miny, maxx, maxy, handles, max_num_handles) iji2dgrid_get_overlaps((self), iji2dgrid_coordinate((self), (minx), (miny)), iji2dgrid_coordinate((self), (maxx), (maxy)), (handles), (max_num_handles))
#define iji2dgrid_get_circle_overlaps(self, centerx, centery, radius, handles, max_num_handles) iji2dgrid_get_overlaps((self), iji2dgrid_coordinate((self), (centerx)-(radius), (centery)-(radius)), iji2dgrid_coordinate((self), (centerx)+(radius), (centery)+(radius)), (handles), (max_num_handles))

/* compute the cell-id from a cell coordinate (x, y) */
#define iji2dgrid_cellid(self, cellx, celly) ((cellx)+(self)->ncols*(celly))
/* extract the column/row (x,y) respectively from the cellid */
#define iji2dgrid_cellid_x(self, cellid) (cellid)%(self)->ncols
#define iji2dgrid_cellid_y(self, cellid) (cellid)/(self)->ncols

/* compute the cell-id coordinates of the cell-id and mask it into a word */
#define iji2dgrid_cellid_to_coordinate(self, cellid) ((iji2dgrid_cellid_x(self, cellid)) | ((iji2dgrid_cellid_y(self, cellid))<<16))

/* returns the x and y of the coordinate respectively */
#define iji2dgrid_coordinate_x(coord) ((coord)&0xffffu)
#define iji2dgrid_coordinate_y(coord) ((coord)>>16)

#define iji2dgrid_coordinate_to_cellid(self, coordinate) iji2dgrid_cellid((self), iji2dgrid_coordinate_x(coordinate), iji2dgrid_coordinate_y(coordinate))

#define iji2dgrid_pointer_add(restype, p, nbytes) (restype)((unsigned char *)(p)+(nbytes))

/* get the data of a handle (or index) (NB: no validity checks) */
#define iji2dgrid_item(self, index_or_handle) iji2dgrid_pointer_add(void*, (self)->data, ((self)->items_mask&(index_or_handle))*(self)->item_size_bytes)

/* get the data of a cellid (NB: no validity checks, i.e. if the grid was initialized with per cell data) */
#define iji2dgrid_cell(restype, self, cellid) iji2dgrid_pointer_add(restype, (self)->celldata, (unsigned)(self)->celldata_size_bytes*(cellid))
#define iji2dgrid_num_cells(self) ((self)->ncols*(self)->nrows)

/*
returns non-zero for success and fills aabb_out with [min_x, min_y, max_x, max_y]
returns zero for failure (invalid cellid), leaving aabb_out untouched
*/
IJI2DGRID_API int iji2dgrid_cell_aabb(struct iji2dgrid *self, unsigned cellid, float aabb_out[4]);

#ifdef __cplusplus
}
#endif

#endif /* IJI2DGRID_INCLUDED_H */

#if defined(IJI2DGRID_IMPLEMENTATION)

#ifndef IJI2DGRID_memset
   #include <string.h>
   #define IJI2DGRID_memset memset
#endif

#ifndef IJI2DGRID_assert
   #include <assert.h>
   #define IJI2DGRID_assert assert
#endif

#ifndef IJI2DGRID_ceilf
   #include <math.h> /* ceil */
   #define IJI2DGRID_ceilf ceilf
#endif

#if _WIN32
   #ifdef __cplusplus
      #define IJI2DG__EXTERNC_DECL_BEGIN extern "C" {
      #define IJI2DG__EXTERNC_DECL_END }
   #else
      #define IJI2DG__EXTERNC_DECL_BEGIN
      #define IJI2DG__EXTERNC_DECL_END
   #endif

   #if defined(_WIN64)
      IJI2DG__EXTERNC_DECL_BEGIN
         unsigned char _BitScanForward64(unsigned long * _Index, unsigned __int64 _Mask);
      IJI2DG__EXTERNC_DECL_END

      typedef unsigned __int64 iji2dgrid_type_t;
      #pragma intrinsic(_BitScanForward64)

      #define iji2dgrid_bitscan_forward _BitScanForward64
      #define IJI2DGRID_TYPE_ALIGMENT 8
   #else
      IJI2DG__EXTERNC_DECL_BEGIN
         unsigned char _BitScanForward(unsigned long * _Index, unsigned long _Mask);
      IJI2DG__EXTERNC_DECL_END

      typedef unsigned long iji2dgrid_type_t;
      #pragma intrinsic(_BitScanForward)

      #define iji2dgrid_bitscan_forward _BitScanForward
      #define IJI2DGRID_TYPE_ALIGMENT 4
   #endif /* _WIN64 */

   static int iji2dgrid__ffs(iji2dgrid_type_t word)
   {
      unsigned long index;
      return iji2dgrid_bitscan_forward(&index, word) ? index : -1;
   }

   #define IJI2DGRID_TYPE_NUM_BITS (sizeof(iji2dgrid_type_t)*8)
   #define IJI2DGRID_TYPE_MASK (IJI2DGRID_TYPE_NUM_BITS-1)

   IJI2DG__EXTERNC_DECL_BEGIN
      unsigned char _BitScanReverse(unsigned long * _Index, unsigned long _Mask);
   IJI2DG__EXTERNC_DECL_END

   #pragma intrinsic(_BitScanReverse)

   static int iji2dgrid__clz32(unsigned val)
   {
      unsigned long idx;
      return _BitScanReverse(&idx, val) ? 31 - idx : 32;
   }
#else
   #error "unsupported platform, ffs and clz implementations missing"
#endif

#define IJI2DGRID_HANDLE_NUM_BITS (sizeof(iji2dgrid_handle)*8)

typedef unsigned int iji2dg_uint32;

#ifdef _MSC_VER
   typedef unsigned __int64 iji2dg_uint64;
#else
   typedef unsigned long long iji2dg_uint64;
#endif

#if defined(__ppc64__) || defined(__aarch64__) || defined(_M_X64) || defined(__x86_64__) || defined(__x86_64)
   typedef iji2dg_uint64 iji2dg_uintptr;
#else
   typedef iji2dg_uint32 iji2dg_uintptr;
#endif

#define iji2dgrid__join2(x, y) x ## y
#define iji2dgrid__join(x, y) iji2dgrid__join2(x, y)
#define iji2dgrid__static_assert(exp) typedef char iji2dgrid__join(static_assert, __LINE__) [(exp) ? 1 : -1]

#if defined(IJI2DGRID_32BIT_HANDLES)
   iji2dgrid__static_assert(sizeof(iji2dgrid_handle) == 4);
#else
   iji2dgrid__static_assert(sizeof(iji2dgrid_handle) == 8);
#endif

#if defined(_WIN64)
   iji2dgrid__static_assert(sizeof(struct iji2dgrid) == 64);
#endif

#define IJI2DGRID_FREELIST_END (unsigned)-1

static iji2dg_uintptr iji2dgrid__nbytes_to_aligned_address(const void *p, unsigned align)
{
   iji2dg_uintptr a = ~(iji2dg_uintptr)p + 1;
   return a & (align - 1);
}

static iji2dgrid_type_t *iji2dgrid__columns(struct iji2dgrid *self)
{
   void *base = iji2dgrid_pointer_add(void*, self->celldata, (self->ncols*self->nrows)*self->celldata_size_bytes);
   return iji2dgrid_pointer_add(iji2dgrid_type_t*, base, iji2dgrid__nbytes_to_aligned_address(base, IJI2DGRID_TYPE_ALIGMENT));
}

static int iji2dgrid__num_bits_per_coord(unsigned ncols, unsigned nrows)
{
   unsigned m = ncols > nrows ? ncols : nrows;
   IJI2DGRID_assert(m);
   m -= 1;
   return m ? (32 - iji2dgrid__clz32(m)) : 0;
}

static int iji2dgrid__clamp(int v, int vmin, int vmax) { return v > vmin ? (v > vmax ? vmax : v) : vmin; }

IJI2DGRID_API unsigned iji2dgrid_coordinate(struct iji2dgrid *self, float posx, float posy)
{
   float gridmin_x = self->gridmin[0], gridmin_y = self->gridmin[1];
   float gridmax_x = self->gridmax[0], gridmax_y = self->gridmax[1];
   float gridextent_x = gridmax_x - gridmin_x;
   float gridextent_y = gridmax_y - gridmin_y;
   float gridlocal_x = posx - gridmin_x;
   float gridlocal_y = posy - gridmin_y;
   int gridcoordx, gridcoordy;

   gridcoordx = (int)((gridlocal_x / gridextent_x) * (float)(self->ncols));
   gridcoordy = (int)((gridlocal_y / gridextent_y) * (float)(self->nrows));
   gridcoordx = iji2dgrid__clamp(gridcoordx, 0, self->ncols-1);
   gridcoordy = iji2dgrid__clamp(gridcoordy, 0, self->nrows-1);
   return gridcoordx | (gridcoordy << 16);
}

IJI2DGRID_API unsigned iji2dgrid_memory_size_needed(float width, float height, float cell_size,
   unsigned max_num_items, unsigned item_size_bytes,
   unsigned cellitem_size_bytes, unsigned cellitem_alignment)
{
   iji2dg_uint64 size_needed = 0;
   iji2dg_uint64 ncols = (iji2dg_uint64)IJI2DGRID_ceilf(width / cell_size);
   iji2dg_uint64 nrows = (iji2dg_uint64)IJI2DGRID_ceilf(height / cell_size);
   iji2dg_uint64 num_grid_words_per_row_col = ((iji2dg_uint64)max_num_items + IJI2DGRID_TYPE_MASK) / IJI2DGRID_TYPE_NUM_BITS;

   item_size_bytes = item_size_bytes < sizeof(unsigned) ? sizeof(unsigned) : item_size_bytes; /* inplace freelist storage requirement */
   IJI2DGRID_assert(iji2dgrid_valid_configuration(width, height, cell_size, max_num_items, item_size_bytes, cellitem_size_bytes, cellitem_alignment));

   size_needed += (ncols*nrows)*cellitem_size_bytes + (cellitem_alignment ? (cellitem_alignment-1) : 0);
   size_needed += ((iji2dg_uint64)max_num_items*item_size_bytes);
   size_needed += (ncols + nrows)*num_grid_words_per_row_col*sizeof(iji2dgrid_type_t);
   size_needed += IJI2DGRID_TYPE_ALIGMENT - 1;

   IJI2DGRID_assert((size_needed&((iji2dg_uint64)-1)<<32)==0);

   return (unsigned)size_needed;
}

IJI2DGRID_API void iji2dgrid_dimensions(float width, float height, float cell_size,
   unsigned *ncols, unsigned *nrows)
{
   unsigned nc = (unsigned)IJI2DGRID_ceilf(width / cell_size);
   unsigned nr = (unsigned)IJI2DGRID_ceilf(height / cell_size);

   IJI2DGRID_assert(nc && nr);

   *ncols = nc;
   *nrows = nr;
}

IJI2DGRID_API unsigned iji2dgrid_max_num_items(float width, float height, float cell_size)
{
   unsigned ncols = (unsigned)IJI2DGRID_ceilf(width / cell_size);
   unsigned nrows = (unsigned)IJI2DGRID_ceilf(height / cell_size);

   int nbits = iji2dgrid__num_bits_per_coord(ncols, nrows);

   nbits *= 4;
   if (nbits >= IJI2DGRID_HANDLE_NUM_BITS)
      return 0;

   nbits = IJI2DGRID_HANDLE_NUM_BITS - nbits;
   if (nbits >= 32)
      return 0xffffffffu;

   return 1u << (32 - nbits);
}

IJI2DGRID_API int iji2dgrid_valid_configuration(float width, float height,
   float cell_size, unsigned max_num_items, unsigned item_size_bytes,
   unsigned cellitem_size_bytes, unsigned cellitem_alignment)
{
   unsigned ncols = (unsigned)IJI2DGRID_ceilf(width / cell_size);
   unsigned nrows = (unsigned)IJI2DGRID_ceilf(height / cell_size);

   int nbits_items = max_num_items > 1 ? (32-iji2dgrid__clz32(max_num_items-1)) : 0;
   int nbits_per_coord = iji2dgrid__num_bits_per_coord(ncols, nrows);

   /* as we encode the index + coordinate_start(x,y) (2 entries) + coordinate_end(x,y) (2 entries) */
   if (0 > ((int)IJI2DGRID_HANDLE_NUM_BITS-nbits_items-4*nbits_per_coord))
      return 0;

   if (2 > max_num_items) /* 2 is minimum size */
      return 0;

   if (0.0f > width)
      return 0;

   if (0.0f > height)
      return 0;

   if (ncols - 1 > 0xffffu) /* we store the columns/rows as unsigned short */
      return 0;

   if (nrows - 1 > 0xffffu) /* we store the columns/rows as unsigned short */
      return 0;

   if (!cellitem_size_bytes != !cellitem_alignment) /* either both zero or both non-zero */
      return 0;

   if (cellitem_size_bytes > 0xffffu) /* we store celldata_size_bytes as unsigned short */
      return 0;

   if (item_size_bytes > 0xffffu) /* we store item_size_bytes as unsigned short */
      return 0;

   return 1;
}

IJI2DGRID_API void iji2dgrid_init(struct iji2dgrid *self,
   float gridminx, float gridminy, float gridmaxx, float gridmaxy,
   float cell_size, unsigned max_num_items, unsigned item_size_bytes,
   unsigned cellitem_size_bytes, unsigned cellitem_alignment, void *memory)
{
   float width = gridmaxx - gridminx;
   float height = gridmaxy - gridminy;
   unsigned ncols = (unsigned)IJI2DGRID_ceilf(width / cell_size);
   unsigned nrows = (unsigned)IJI2DGRID_ceilf(height / cell_size);
   int nbits_items = max_num_items > 1 ? (32-iji2dgrid__clz32(max_num_items-1)) : 0;

   IJI2DGRID_assert(iji2dgrid_valid_configuration(width, height, cell_size, max_num_items, item_size_bytes, cellitem_size_bytes, cellitem_alignment));

   self->ncols = (unsigned short)ncols;
   self->nrows = (unsigned short)nrows;
   self->nbits_per_coord = (unsigned char)iji2dgrid__num_bits_per_coord(ncols, nrows);
   self->nbits_items = (unsigned char)nbits_items;

   self->cell_size = cell_size;
   self->gridmin[0] = gridminx;
   self->gridmin[1] = gridminy;
   self->gridmax[0] = gridmaxx;
   self->gridmax[1] = gridmaxy;

   self->items_capacity = max_num_items;
   self->items_mask = 0xffffffffu >> (32 - nbits_items);
   self->item_size_bytes = (unsigned short)item_size_bytes;

   self->ngridwords_per_row_col = ((iji2dg_uint64)max_num_items + IJI2DGRID_TYPE_MASK) / IJI2DGRID_TYPE_NUM_BITS;

   self->data = memory;
   self->celldata = iji2dgrid_pointer_add(void*, self->data, self->items_capacity*self->item_size_bytes);
   if (cellitem_alignment)
      self->celldata = iji2dgrid_pointer_add(void*, self->celldata, iji2dgrid__nbytes_to_aligned_address(self->celldata, cellitem_alignment));

   self->celldata_size_bytes = (unsigned short)cellitem_size_bytes;
   self->feature_flags = 0;

   if (32 >= self->nbits_items + 4 * self->nbits_per_coord)
      self->feature_flags |= IJI2DGRID_FF_32BIT_CLEAN_HANDLES;

   iji2dgrid_reset(self);
}

IJI2DGRID_API void iji2dgrid_reset(struct iji2dgrid *self)
{
   unsigned i, items_capacity = self->items_capacity;
   unsigned item_size_bytes = self->item_size_bytes;
   unsigned total_num_gridwords;
   void *p = self->data;

   self->data_freelist = 1; /* so we can ensure that 0 is a invalid handle */

   for (i = 0; i != items_capacity - 1; ++i) {
      *(unsigned*)p = i + 1;
      p = iji2dgrid_pointer_add(void*, p, item_size_bytes);
   }
   *(unsigned*)p = IJI2DGRID_FREELIST_END;

   p = (void*)iji2dgrid__columns(self);

   total_num_gridwords = (self->ncols + self->nrows) * self->ngridwords_per_row_col;
   IJI2DGRID_memset(p, 0, sizeof(iji2dgrid_type_t)*total_num_gridwords);
}

#define iji2dgrid__encode_coord(packhandle_nbits, nbits_per_coord, x, y, xend, yend) (((iji2dgrid_handle)(x)<<(3*nbits_per_coord))|((iji2dgrid_handle)(y)<<(2*nbits_per_coord))|((iji2dgrid_handle)(xend)<<(1*nbits_per_coord))|((iji2dgrid_handle)(yend)<<(0*nbits_per_coord)))<<((packhandle_nbits)-4*nbits_per_coord)

static void iji2dgrid__decompose_handle(struct iji2dgrid *self, iji2dgrid_handle handle, unsigned *x, unsigned *y, unsigned *xend, unsigned *yend)
{
   unsigned nbits_per_coord = self->nbits_per_coord;
   unsigned packhandle_nbits = (self->feature_flags&IJI2DGRID_FF_32BIT_CLEAN_HANDLES) ? 32 : IJI2DGRID_HANDLE_NUM_BITS;
   iji2dgrid_handle coordinate_mask = (iji2dgrid_handle)-1 << (IJI2DGRID_HANDLE_NUM_BITS - nbits_per_coord);
   unsigned shift = packhandle_nbits - nbits_per_coord;

   if (self->feature_flags&IJI2DGRID_FF_32BIT_CLEAN_HANDLES)
      coordinate_mask >>= (IJI2DGRID_HANDLE_NUM_BITS-32);

   *x = (unsigned)((handle & coordinate_mask) >> shift);
   coordinate_mask >>= nbits_per_coord;
   shift -= nbits_per_coord;

   *y = (unsigned)((handle & coordinate_mask) >> shift);
   coordinate_mask >>= nbits_per_coord;
   shift -= nbits_per_coord;

   *xend = (unsigned)((handle & coordinate_mask) >> shift);
   coordinate_mask >>= nbits_per_coord;
   shift -= nbits_per_coord;

   *yend = (unsigned)((handle & coordinate_mask) >> shift);
}

IJI2DGRID_API void *iji2dgrid_add_box(struct iji2dgrid *self, float minx, float miny, float maxx, float maxy, iji2dgrid_handle *handle_out)
{
   unsigned free_list = self->data_freelist;
   unsigned coord_start, coord_end;
   unsigned x, xstart, xend, y, ystart, yend;
   unsigned index_bucket, index_position;
   void *res;
   iji2dgrid_type_t *gc, *gr;
   unsigned item_index, ngridwords_per_row_col = self->ngridwords_per_row_col, nbits_per_coord = self->nbits_per_coord;
   unsigned packhandle_nbits = (self->feature_flags&IJI2DGRID_FF_32BIT_CLEAN_HANDLES) ? 32 : IJI2DGRID_HANDLE_NUM_BITS;

   if (free_list == IJI2DGRID_FREELIST_END) {
      *handle_out = 0;
      return 0;
   }

   IJI2DGRID_assert(maxx >= minx && maxy >= miny);

   IJI2DGRID_assert((free_list & ~self->items_mask) == 0);
   res = iji2dgrid_pointer_add(void *, self->data, free_list*self->item_size_bytes);
   self->data_freelist = *(unsigned*)res;

   item_index = free_list;
   IJI2DGRID_assert((item_index & ~self->items_mask) == 0);

   index_bucket = item_index / IJI2DGRID_TYPE_NUM_BITS;
   index_position = item_index & IJI2DGRID_TYPE_MASK;

   coord_start = iji2dgrid_coordinate(self, minx, miny);
   coord_end = iji2dgrid_coordinate(self, maxx, maxy);

   gc = iji2dgrid__columns(self);
   x = xstart = iji2dgrid_coordinate_x(coord_start);
   for (xend = iji2dgrid_coordinate_x(coord_end); x <= xend; ++x) {
      iji2dgrid_type_t *c = gc + x*ngridwords_per_row_col;

      c[index_bucket] |= (iji2dgrid_type_t)1 << index_position;
   }

   gr = gc + self->ncols*self->ngridwords_per_row_col;
   y = ystart = iji2dgrid_coordinate_y(coord_start);
   for (yend = iji2dgrid_coordinate_y(coord_end); y <= yend; ++y) {
      iji2dgrid_type_t *r = gr + y*ngridwords_per_row_col;

      r[index_bucket] |= (iji2dgrid_type_t)1 << index_position;
   }

   *handle_out = (iji2dgrid_handle)item_index | (iji2dgrid_handle)iji2dgrid__encode_coord(packhandle_nbits, nbits_per_coord, xstart, ystart, xend, yend);
   #if !defined(IJI2DGRID_32BIT_HANDLES)
      IJI2DGRID_assert(!(self->feature_flags&IJI2DGRID_FF_32BIT_CLEAN_HANDLES) || ((*handle_out) >> 32) == 0);
   #endif

   return res;
}

/*
remove item from grid where the coordinates is not encoded in the handle.
as we currently does not support _not_ encoding grid coordinates in the handles
let us not expose this in the public api.
*/
static void iji2dgrid__removeex(struct iji2dgrid *self, iji2dgrid_handle handle,
   unsigned x, unsigned y, unsigned xend, unsigned yend)
{
   unsigned ngridwords_per_row_col = self->ngridwords_per_row_col;
   unsigned  item_index;
   unsigned index_bucket, index_position;
   iji2dgrid_type_t *gc, *gr;
   void *p;

   if (!handle)
      return;

   IJI2DGRID_assert(self->ncols > x && self->ncols > xend && xend >= x);
   IJI2DGRID_assert(self->nrows > y && self->nrows > yend && yend >= y);

   item_index = handle&self->items_mask;
   IJI2DGRID_assert(self->items_capacity > item_index); /* can happen as we do not force power of two num items */
   index_bucket = item_index / IJI2DGRID_TYPE_NUM_BITS;
   index_position = item_index & IJI2DGRID_TYPE_MASK;

   p = iji2dgrid_item(self, item_index);
   *(unsigned*)p = self->data_freelist;
   self->data_freelist = item_index;

   gc = iji2dgrid__columns(self);
   for (; x <= xend; ++x) {
      iji2dgrid_type_t *c = gc + x*ngridwords_per_row_col;

      c[index_bucket] &= ~((iji2dgrid_type_t)1 << index_position);
   }

   gr = gc + self->ncols*ngridwords_per_row_col;
   for (; y <= yend; ++y) {
      iji2dgrid_type_t *r = gr + y*ngridwords_per_row_col;

      r[index_bucket] &= ~((iji2dgrid_type_t)1 << index_position);
   }
}

IJI2DGRID_API void iji2dgrid_remove(struct iji2dgrid *self, iji2dgrid_handle handle)
{
   unsigned x, y, xend, yend;
   iji2dgrid__decompose_handle(self, handle, &x, &y, &xend, &yend);
   iji2dgrid__removeex(self, handle, x, y, xend, yend);
}

IJI2DGRID_API unsigned iji2dgrid_get_overlaps(struct iji2dgrid *self, unsigned coord_start, unsigned coord_end, void **handles, unsigned max_num_handles)
{
   unsigned outcount = 0;
   iji2dgrid_type_t *gc, *gr;
   unsigned ngridwords_per_row_col = self->ngridwords_per_row_col;

   gc = iji2dgrid__columns(self);
   gr = gc + self->ncols*ngridwords_per_row_col;

   if (coord_start == coord_end) {
      unsigned i;
      iji2dgrid_type_t *c = gc + iji2dgrid_coordinate_x(coord_start)*ngridwords_per_row_col;
      iji2dgrid_type_t *r = gr + iji2dgrid_coordinate_y(coord_start)*ngridwords_per_row_col;
      IJI2DGRID_assert(iji2dgrid_coordinate_x(coord_start) < self->ncols);
      IJI2DGRID_assert(iji2dgrid_coordinate_y(coord_start) < self->nrows);

      for (i = 0; i != ngridwords_per_row_col; ++i) {
         iji2dgrid_type_t all_items_mask = c[i] & r[i];

         while (all_items_mask) {
            /*
            * suppress MSVC's warning 4146 (unary minus operator applied to unsigned type, result still unsigned)
            * by doing it 'manually'.
            * should read:
            *     iji2dgrid_type_t item_mask = all_items_mask & -all_items_mask;
            */
            iji2dgrid_type_t item_mask = all_items_mask & ((iji2dgrid_type_t)0-all_items_mask); /* isolate rightmost bit */

            if (max_num_handles != 0) {
               int index = iji2dgrid__ffs(item_mask);

               handles[outcount] = iji2dgrid_pointer_add(void*, self->data, (i*IJI2DGRID_TYPE_NUM_BITS + index)*self->item_size_bytes);
               --max_num_handles;
            }

            ++outcount;

            all_items_mask ^= item_mask;
         }
      }
   } else {
      unsigned i, x, y;
      unsigned mincx, mincy, maxcx, maxcy;
      mincx = maxcx = iji2dgrid_coordinate_x(coord_start);
      mincy = maxcy = iji2dgrid_coordinate_y(coord_start);

      mincx = mincx > iji2dgrid_coordinate_x(coord_end) ? iji2dgrid_coordinate_x(coord_end) : mincx;
      maxcx = maxcx < iji2dgrid_coordinate_x(coord_end) ? iji2dgrid_coordinate_x(coord_end) : maxcx;

      mincy = mincy > iji2dgrid_coordinate_y(coord_end) ? iji2dgrid_coordinate_y(coord_end) : mincy;
      maxcy = maxcy < iji2dgrid_coordinate_y(coord_end) ? iji2dgrid_coordinate_y(coord_end) : maxcy;

      IJI2DGRID_assert(mincx < self->ncols && maxcx < self->ncols);
      IJI2DGRID_assert(mincy < self->nrows && maxcy < self->nrows);

      for (i = 0; i != ngridwords_per_row_col; ++i) {
         iji2dgrid_type_t all_items_mask, all_items_mask_col = 0, all_items_mask_row = 0;

         for (x = mincx; x <= maxcx; ++x) {
            iji2dgrid_type_t *c = gc + x*ngridwords_per_row_col;
            all_items_mask_col |= c[i];
         }

         for (y = mincy; y <= maxcy; ++y) {
            iji2dgrid_type_t *r = gr + y*ngridwords_per_row_col;

            all_items_mask_row |= r[i];
         }

         all_items_mask = all_items_mask_col & all_items_mask_row;

         while (all_items_mask) {
            iji2dgrid_type_t item_mask = all_items_mask & ((iji2dgrid_type_t)0-all_items_mask); /* isolate rightmost bit */

            if (max_num_handles != 0) {
               int index = iji2dgrid__ffs(item_mask);

               handles[outcount] = iji2dgrid_pointer_add(void*, self->data, (i*IJI2DGRID_TYPE_NUM_BITS + index)*self->item_size_bytes);
               --max_num_handles;
            }

            ++outcount;

            all_items_mask ^= item_mask;
         }
      }
   }

   return outcount;
}

IJI2DGRID_API int iji2dgrid_cell_aabb(struct iji2dgrid *self, unsigned cellid, float aabb_out[4])
{
   unsigned cellx, celly;

   if (cellid >= (unsigned)self->ncols*self->nrows)
      return 0;

   /* as we reliable only can check celly, as we use % for cellx, we assume this is a valid cellid */
   cellx = iji2dgrid_cellid_x(self, cellid);
   celly = iji2dgrid_cellid_y(self, cellid);

   aabb_out[0] = self->gridmin[0] + (float)cellx*self->cell_size;
   aabb_out[1] = self->gridmin[1] + (float)celly*self->cell_size;

   aabb_out[2] = aabb_out[0] + self->cell_size;
   aabb_out[3] = aabb_out[1] + self->cell_size;

   aabb_out[2] = (aabb_out[2] > self->gridmax[0] ? self->gridmax[0] : aabb_out[2]);
   aabb_out[3] = (aabb_out[3] > self->gridmax[1] ? self->gridmax[1] : aabb_out[3]);

   return 1;
}

#if defined (IJI2DGRID_TEST) || defined (IJI2DGRID_TEST_MAIN)

char iji2dgrid__gridstorage[1024*16];

struct iji2dgrid_userdata {
   iji2dgrid_handle handle;
   unsigned i;
   unsigned coord;
};

static void iji2dgrid_test(void)
{
#define IJ2DGRID_MAX_NUM_OBJECTS (80)
   struct iji2dgrid sgrid, *self = &sgrid;
   float gridextents[2] = { 2047.5f, 2047.5f };
   float cell_size = 128.0f;
   unsigned c, i, item_size = sizeof(struct iji2dgrid_userdata);
   unsigned max_num_items = IJ2DGRID_MAX_NUM_OBJECTS;
   iji2dgrid_handle all_handles[IJ2DGRID_MAX_NUM_OBJECTS];

   void *ps[IJ2DGRID_MAX_NUM_OBJECTS];

   int loop;
   unsigned origocoord, origocx, origocy;
   float width, height, aabb[4];
   unsigned x, y;
   aabb[0] = -gridextents[0] * 0.5f;
   aabb[1] = -gridextents[0] * 0.5f;

   aabb[2] = aabb[0] + gridextents[0];
   aabb[3] = aabb[1] + gridextents[1];
   width = aabb[2] - aabb[0];
   height = aabb[3] - aabb[1];
   IJI2DGRID_assert((sizeof ps / sizeof *ps) >= max_num_items);

   for (loop = 0; loop != 4; ++loop) {
      unsigned s, stores, max_stores = 3;
      float aabb_origo[4];
      unsigned cellitem_size = loop*sizeof(unsigned);
      unsigned cellitem_alignment = loop ? 4 : 0;
      IJI2DGRID_memset(all_handles, 0, sizeof all_handles);
      item_size += loop;
      max_num_items += loop*IJI2DGRID_TYPE_NUM_BITS;
      IJI2DGRID_assert(iji2dgrid_memory_size_needed(width, height, cell_size, max_num_items, item_size, cellitem_size, cellitem_alignment) <= sizeof(iji2dgrid__gridstorage));
      iji2dgrid_init(self, aabb[0], aabb[1], aabb[2], aabb[3], cell_size, max_num_items, item_size,
          cellitem_size, cellitem_alignment,
          (void*)iji2dgrid__gridstorage);

      origocoord = iji2dgrid_coordinate(self, 0.0f, 0.0f);
      origocx = origocoord & 0xffffu;
      origocy = origocoord >> 16;

      if (cellitem_size >= 4) {
         for (x = 0; x != self->ncols; ++x) {
            for (y = 0; y != self->nrows; ++y) {
               unsigned coord = iji2dgrid_cellid(self, x, y);
               unsigned *p = iji2dgrid_cell(unsigned*, self, coord);
               *p = coord;
            }
         }
         /* now verify our stores above */
         for (x = 0; x != self->ncols; ++x) {
            for (y = 0; y != self->nrows; ++y) {
               unsigned coord = iji2dgrid_cellid(self, x, y);
               unsigned *p = iji2dgrid_cell(unsigned*, self, coord);

               IJI2DGRID_assert(*(unsigned*)p == coord);
            }
         }
      }

       /* make make an overlap over whole grid */
      c = iji2dgrid_get_overlaps(self, 0, ((self->ncols-1) | ((unsigned)(self->nrows-1))<<16), ps, sizeof ps / sizeof *ps);
      IJI2DGRID_assert(c == 0);

      iji2dgrid_cell_aabb(self, iji2dgrid_coordinate_to_cellid(self, origocoord), aabb_origo);

      for (stores = 0; stores != max_stores; ++stores) {
         iji2dgrid_handle temp_handle_storage[4];
         unsigned num_items_to_add = stores+1;
         iji2dgrid_handle handle, origohandle;
         struct iji2dgrid_userdata *r;
         unsigned this_coord = origocoord, this_y = origocy;
         unsigned minycoord = origocy, maxycoord = origocy;
         unsigned nactive;
         float cx, cy, cr = 1.0f;
         unsigned prev_coord = origocoord;
         IJI2DGRID_assert((sizeof temp_handle_storage/sizeof *temp_handle_storage) > num_items_to_add);
         cx = (aabb_origo[0] + aabb_origo[2])*0.5f;
         cy = (aabb_origo[1] + aabb_origo[3])*0.5f;
         r = (struct iji2dgrid_userdata*)iji2dgrid_add_cirle(self, cx, cy, cr, &handle);
         IJI2DGRID_assert(r != 0);
         origohandle = r->handle = handle;
         r->coord = origocoord;
         c = iji2dgrid_get_overlaps(self, origocoord, origocoord, ps, sizeof ps / sizeof *ps);
         IJI2DGRID_assert(c == 1);
         IJI2DGRID_assert(((struct iji2dgrid_userdata*)ps[0])->handle == origohandle);

         nactive = 1;

         /* now add 'stores' num item(s) to each row 'above' origo */
         while (this_y != self->nrows-1u) {
            unsigned check, mincoord, maxcoord;
            cy += self->cell_size;
            this_coord = iji2dgrid_coordinate(self, cx, cy);
            IJI2DGRID_assert(this_coord != prev_coord);
            prev_coord = this_coord;
            this_y = this_coord>>16;
            minycoord = minycoord > this_y ? this_y : minycoord;
            maxycoord = maxycoord > this_y ? maxycoord : this_y;
            for (s = 0; s != num_items_to_add; ++s) {
               r = (struct iji2dgrid_userdata*)iji2dgrid_add_cirle(self, cx, cy, cr, &handle);
               IJI2DGRID_assert(r != 0);
               r->handle = handle;
               r->coord = this_coord;
               temp_handle_storage[s] = handle;
               /* if you want auxiliary data with you handles you can
               get the sparse index by or:ing with the grids items-mask */
               all_handles[handle&self->items_mask] = handle;
               ++nactive;
            }

            c = iji2dgrid_get_overlaps(self, this_coord, this_coord, ps, sizeof ps / sizeof *ps);
            IJI2DGRID_assert(c == num_items_to_add);
            for (s = 0; s != num_items_to_add; ++s) {
               iji2dgrid_handle handle_to_find = ((struct iji2dgrid_userdata*)ps[s])->handle;
               int found = 0;
               for (check = 0; check != num_items_to_add; ++check) {
                  if (temp_handle_storage[check] == handle_to_find) {
                     temp_handle_storage[check] = 0;
                     found = 1;
                  }
               }
               IJI2DGRID_assert(found);
            }

            c = iji2dgrid_get_overlaps(self, origocoord, origocoord, ps, sizeof ps / sizeof *ps);
            IJI2DGRID_assert(c == 1);
            IJI2DGRID_assert(((struct iji2dgrid_userdata*)ps[0])->handle == origohandle);
            mincoord = origocx | (minycoord<<16);
            maxcoord = origocx | (maxycoord<<16);

            c = iji2dgrid_get_overlaps(self, mincoord, maxcoord, ps, sizeof ps / sizeof *ps);
            IJI2DGRID_assert(c == nactive);
         }

         c = iji2dgrid_get_box_overlaps(self, self->gridmin[0], self->gridmin[1], self->gridmax[0], self->gridmax[1], ps, sizeof ps / sizeof *ps);
         IJI2DGRID_assert(c == nactive);

         c = iji2dgrid_get_overlaps(self, origocoord, prev_coord, ps, sizeof ps / sizeof *ps);
         IJI2DGRID_assert(c == nactive);

         /* go over all columns and verify overlaps, doing an overlap for each column of all its rows */
         for (i = 0; i != self->ncols; ++i) {
            unsigned mincoord, maxcoord;
            unsigned nexpected = 0;
            if (i == origocx)
               nexpected = nactive;

            mincoord = i;
            maxcoord = i | ((unsigned)self->nrows-1)<<16;
            c = iji2dgrid_get_overlaps(self, mincoord, maxcoord, ps, sizeof ps / sizeof *ps);
            IJI2DGRID_assert(c == nexpected);
         }

         /* go over all rows and verify overlaps, doing an overlap for each row of all its columns */
         for (i = 0; i != self->nrows; ++i) {
            unsigned mincoord, maxcoord;
            unsigned nexpected = 0;
            this_coord = origocx | i<<16;
            if (i == origocy)
               nexpected = 1;
            else if (i > origocy)
               nexpected = num_items_to_add;

            mincoord = 0 | (i<<16);
            maxcoord = (self->ncols-1) | (i<<16);
            c = iji2dgrid_get_overlaps(self, mincoord, maxcoord, ps, sizeof ps / sizeof *ps);
            IJI2DGRID_assert(c == nexpected);
            IJI2DGRID_assert(!c || ((struct iji2dgrid_userdata*)ps[0])->coord == this_coord);
         }

         this_coord = (self->ncols-1) | ((unsigned)(self->nrows-1))<<16;
         c = iji2dgrid_get_overlaps(self, 0, this_coord, ps, sizeof ps / sizeof *ps);
         IJI2DGRID_assert(c == nactive);

         for (i = 0; i != c; ++i) {
            struct iji2dgrid_userdata *fetched = (struct iji2dgrid_userdata *)ps[i];

            iji2dgrid_remove(self, fetched->handle);
         }

         c = iji2dgrid_get_overlaps(self, 0, this_coord, ps, sizeof ps / sizeof *ps);
         IJI2DGRID_assert(c == 0);
      }

      /* now verify our initial stores is unchanged */
      if (cellitem_size >= 4) {
         for (x = 0; x != self->ncols; ++x) {
            for (y = 0; y != self->nrows; ++y) {
               unsigned coord = iji2dgrid_cellid(self, x, y);
               IJI2DGRID_assert(*iji2dgrid_cell(unsigned*, self, coord) == coord);
            }
         }
      }

   }

#undef IJ2DGRID_MAX_NUM_OBJECTS
}

#if defined(IJI2DGRID_TEST_MAIN)
int main(int args, char **argc)
{
   iji2dgrid_test();
   (void)args;
   (void)argc;
   return 0;
}
#endif

#endif /* defined (IJI2DGRID_TEST) || defined (IJI2DGRID_TEST_MAIN) */

#endif /* IJI2DGRID_IMPLEMENTATION */

/*
LICENSE
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - 3-Clause BSD License
Copyright (c) 2018-, Fredrik Engkvist
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
* Neither the name of the copyright holder nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
*/
/* clang-format on */
