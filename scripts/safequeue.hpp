template <class T>
class SafeQueue{
    private:
        std::queue<T> q;
        std::mutex m;
    public:
        void push(T t){
            std::lock_guard<std::mutex> lock(m);
            q.push(t);
        }
        T pop(){
            std::lock_guard<std::mutex> lock(m);
            if ( q.empty() )
                return T();
            T val = q.front();
            q.pop();
            return val;
        }
        int size(){
            std::lock_guard<std::mutex> lock(m);
            return q.size();
        }
        bool empty(){
            std::lock_guard<std::mutex> lock(m);
            return q.empty();
        }
};