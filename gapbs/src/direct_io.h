#ifndef DIRECT_IO_H_
#define DIRECT_IO_H_

#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <fcntl.h>
#include <unistd.h>

class DirectFile {
 public:
  static constexpr size_t kAlign   = 4096;
  static constexpr size_t kBufSize = 1 << 20;  // 1 MiB

  explicit DirectFile(const std::string &filename) : fd_(-1), buf_(nullptr),
      buf_start_(0), buf_end_(0), file_offset_(0) {
    fd_ = open(filename.c_str(), O_RDONLY | O_DIRECT);
    if (fd_ < 0) {
      std::cout << "DirectFile: cannot open " << filename << ": "
                << strerror(errno) << std::endl;
      std::exit(-2);
    }
    if (posix_memalign(reinterpret_cast<void**>(&buf_), kAlign, kBufSize) != 0) {
      std::cout << "DirectFile: posix_memalign failed" << std::endl;
      std::exit(-2);
    }
  }

  ~DirectFile() {
    if (fd_ >= 0) close(fd_);
    free(buf_);
  }

  bool is_open() const { return fd_ >= 0; }

  void read_exact(void *dst, size_t nbytes) {
    uint8_t *out = reinterpret_cast<uint8_t*>(dst);
    while (nbytes > 0) {
      if (buf_start_ == buf_end_)
        refill();
      size_t avail = buf_end_ - buf_start_;
      size_t take  = nbytes < avail ? nbytes : avail;
      std::memcpy(out, buf_ + buf_start_, take);
      buf_start_ += take;
      out        += take;
      nbytes     -= take;
    }
  }

  template<typename T>
  void read_val(T &val) {
    read_exact(&val, sizeof(T));
  }

 private:
  void refill() {
    size_t rounded = kBufSize;
    ssize_t got = ::read(fd_, buf_, rounded);
    if (got < 0) {
      std::cout << "DirectFile: read error: " << strerror(errno) << std::endl;
      std::exit(-2);
    }
    file_offset_ += static_cast<size_t>(got);
    buf_start_ = 0;
    buf_end_   = static_cast<size_t>(got);
  }

  int      fd_;
  uint8_t *buf_;
  size_t   buf_start_;
  size_t   buf_end_;
  size_t   file_offset_;
};

#endif  // DIRECT_IO_H_
