
namespace Mantid {
namespace DataObjects {
template <typename T>
class EventListBaseFunctionsTemplate
{
  

private:
    friend T;
    EventListTemplate() = default;

    inline T & as_underlying()
    {
        return static_cast<T&>(*this);
    }
};
}
}