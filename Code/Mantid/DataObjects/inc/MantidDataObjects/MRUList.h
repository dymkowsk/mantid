#ifndef MRULIST_H
#define MRULIST_H
/*
 * MRUList.h
 *
 *  Created on: Jul 14, 2010
 *      Author: janik
 */

#include "MantidKernel/System.h"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

#include <fstream>
#include <valarray>


namespace Mantid
{
namespace DataObjects
{

  using namespace boost::multi_index;

  /** An MRU (most recently used) list keeps record of the last n
  *  inserted items, listing first the newer ones. Care has to be
  *  taken when a duplicate item is inserted: instead of letting it
  *  appear twice, the MRU list relocates it to the first position.
  *  This class has been taken from one of the examples given in the
  *  Boost.MultiIndex documentation (<http://www.boost.org/libs/multi_index/doc/reference/index.html>)
  */
  template<class T>
  class DLLExport MRUList
  {
  public:

  private:
    /// hideous typedef for the container holding the list
    typedef typename boost::multi_index::multi_index_container<
      T*,
      boost::multi_index::indexed_by<
        boost::multi_index::sequenced<>,
        boost::multi_index::hashed_unique<BOOST_MULTI_INDEX_CONST_MEM_FUN(T,int,minIndex)>
      >
    > item_list;

    /// This typedef makes an ordered item list (you access it by the 1st index)
    typedef typename nth_index<item_list,1>::type ordered_item_list;

    /// The most recently used list
    mutable item_list il;
    /// The length of the list
    const std::size_t max_num_items;

  public:
    /** Constructor
    *  @param max_num_items_ The length of the list
    */
    MRUList(const std::size_t &max_num_items_)  :
      max_num_items(max_num_items_)
      {
      }

  public:
    /** Insert an item into the list. If it's already in the list, it's moved to the top.
    *  If it's a new item, it's put at the top and the last item in the list is written to file and dropped.
    *  @param item The item, of type T, to put in the list
    *  @return pointer to an item that is being dropped from the MRU. The calling code can
    *     do stuff to it (save it) and needs to delete it. NULL if nothing needs to be dropped.
    */
    T* insert(T* item)
    {
      std::pair<typename MRUList<T>::item_list::iterator,bool> p=this->il.push_front(item);

      if (!p.second)
      { /* duplicate item */
        this->il.relocate(this->il.begin(), p.first); /* put in front */
        return NULL;
      }
      else if (this->il.size()>max_num_items)
      { /* keep the length <= max_num_items */
        // This is dropping an item - need to write it to disk (if it's changed) and delete
        // but this is left up to the calling class to do.
        T *toWrite = this->il.back();
        this->il.pop_back();
        return toWrite;
      }
      return NULL;
    }


    /// Delete all the data blocks pointed to by the list, and empty the list itself
    void clear()
    {
      for (typename MRUList<T>::item_list::iterator it = this->il.begin(); it != this->il.end(); ++it)
      {
        delete *it;
      }
      this->il.clear();
    }

    /// Size of the list
    size_t size() const {return il.size();}


  public:
    /** Find an element of the list from the key of the minIndex
    *  @param minIndex The minIndex value to search the list for
    */
    T* find(const unsigned int minIndex) const
    {
      typename ordered_item_list::const_iterator it = il.get<1>().find(minIndex);
      if (it == il.get<1>().end())
      {
        return NULL;
      }
      else
      {
        return *it;
      }
    }

  };



}
}


#endif // MRULIST_H
