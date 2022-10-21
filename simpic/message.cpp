// Message class
//
// Original by HyunChul Kim, 2005

#include "ostring.hpp"
#include "message.hpp"
#include <stdarg.h>

Message::Message()
{
  msg_title = "";
  stream = stderr;
  global_ml = MSG_MED;
  file_opened = false;
}

Message::~Message()
{
  if (file_opened) fclose(stream);
}

void Message::set_msg_level_high()
{
  global_ml = MSG_HIGH;
}

void Message::set_msg_level_med()
{
  global_ml = MSG_MED;
}

void Message::set_msg_level_low()
{
  global_ml = MSG_LOW;
}

void Message::set_msg(const char *_msg_title,
		      const char *filename, Msg_Level _global_ml)
{
  msg_title = _msg_title;
  global_ml = _global_ml;
  stream = fopen(filename, "w");
  if (stream) file_opened = true;
  else
  {
    fprintf(stderr, "\nThe file `%s' cannot be opened.", filename);
    file_opened = false;
  }
}

void Message::set_msg(const char *_msg_title, Msg_Level _global_ml,
		      FILE *_stream)
{
  msg_title = _msg_title;
  global_ml = _global_ml;
  stream = _stream;
}

void Message::error_msg(const char *format, ...) const
{
  va_list argptr;

  print_title();
  va_start(argptr, format);
  vfprintf(stream, format, argptr);
  va_end(argptr);
}

// With title
void Message::print_msg(const char *format, ...) const
{
  if (global_ml >= MSG_MED)
  {
    print_title();
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stream, format, argptr);
    va_end(argptr);
  }
}

// With title
void Message::print_msg(Msg_Level _local_ml, const char *format, ...) const
{
  if (_local_ml+global_ml >= MSG_MED)
  {
    print_title();
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stream, format, argptr);
    va_end(argptr);
  }
}

// Without title
void Message::print_shr_msg(const char *format, ...) const
{
  if (global_ml >= MSG_MED)
  {
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stream, format, argptr);
    va_end(argptr);
  }
}

// Without title
void Message::print_shr_msg(Msg_Level _local_ml, const char *format, ...) const
{
  if (_local_ml+global_ml >= MSG_MED)
  {
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stream, format, argptr);
    va_end(argptr);
  }
}

int Message::terminate_run(const char *format, ...) const
{
  va_list argptr;

  print_title();
  va_start(argptr, format);
  vfprintf(stream, format, argptr);
  va_end(argptr);
  print_shr_msg(MSG_MED,"\nexiting...\n");
  exit(1);
}

