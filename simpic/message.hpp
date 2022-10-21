// Message class
//
// Original by HyunChul Kim, 2005

#ifndef MESSAGE_H
#define MESSAGE_H
#include "ostring.hpp"
#include <stdlib.h>

#ifndef Msg_Level
#ifndef MSG_LEVEL_DECLARED
enum Msg_Level {MSG_LOW=-1,MSG_MED=0,MSG_HIGH=1};
#define MSG_LEVEL_DECLARED
#endif
#endif

class Message
{
private:
  ostring msg_title;
  FILE *stream;
  Msg_Level global_ml;
  bool file_opened;

protected:
  Message();
  ~Message();

  // set how much the given message is important (say, priority).
  void set_msg_level_high();
  void set_msg_level_low();
  void set_msg_level_med();

  // set msg_title to be printed at the front of message.
  void set_msg(const char *_msg_title, Msg_Level _global_ml=MSG_MED,
	       FILE *_stream=stderr);
  void set_msg(const char *_msg_title, const char *filename,
	       Msg_Level _global_ml=MSG_MED);

  // terminate_run
  int terminate_run(const char *format, ...) const;
  // print out error message but not die
  void error_msg(const char *format, ...) const;

  inline void print_title() const
  {
    fprintf(stream, "\n%s:: ", msg_title());
  }

  // print out the message (information rather than error)
  void print_msg(const char *format, ...) const;
  void print_shr_msg(const char *format, ...) const;
  void print_msg(Msg_Level _local_ml, const char *format, ...) const;
  void print_shr_msg(Msg_Level _local_ml, const char *format, ...) const;
};
#endif
