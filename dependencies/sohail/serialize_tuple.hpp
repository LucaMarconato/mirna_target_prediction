// origin: http://uint32t.blogspot.com/2008/03/serializing-boosttuple-using.html

#include <boost/tuple/tuple.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/preprocessor/repetition.hpp>

namespace boost { namespace serialization {

    template<typename Archive, typename T1>

    void serialize(Archive & ar,
                   boost::tuples::cons<T1,boost::tuples::null_type> & t,

                   const unsigned int)
    {
      ar & boost::serialization::make_nvp("head",t.head);
    }

    template<typename Archive, typename T1, typename T2>

    void serialize(Archive & ar,
                   boost::tuples::cons<T1,T2> & t,

                   const unsigned int)
    {
      ar & boost::serialization::make_nvp("head",t.head);

      ar & boost::serialization::make_nvp("tail",t.tail);
    }

    template<typename Archive, typename T1>
    void serialize(Archive & ar,

                   boost::tuple<T1> & t,
                   const unsigned int)
    {
      ar & boost::serialization::make_nvp("head",t.head);
    }

#define GENERATE_TUPLE_SERIALIZE(z,nargs,unused)                            \
    template< typename Archive, BOOST_PP_ENUM_PARAMS(nargs,typename T) > \
    void serialize(Archive & ar,                                        \
                   boost::tuple< BOOST_PP_ENUM_PARAMS(nargs,T) > & t,   \
                   const unsigned int)                          \
    {                                                                   \
      ar & boost::serialization::make_nvp("head",t.head);               \
      ar & boost::serialization::make_nvp("tail",t.tail);               \
    }


    BOOST_PP_REPEAT_FROM_TO(2,6,GENERATE_TUPLE_SERIALIZE,~);

}}
