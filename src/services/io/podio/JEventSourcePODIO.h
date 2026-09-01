// Copyright 2022, David Lawrence
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <cstdint>
#include <JANA/JApplicationFwd.h>
#include <JANA/JEventSource.h>
#include <JANA/JEventSourceGeneratorT.h>
#include <podio/Frame.h>
#include <podio/Reader.h>
#include <spdlog/logger.h>
#include <cstddef>
#include <memory>
#include <chrono>
#include <set>
#include <string>
#include <string_view>
#include <vector>

class JEventSourcePODIO : public JEventSource {

public:
  JEventSourcePODIO();

  virtual ~JEventSourcePODIO();

  void Init() override;

  void Open() override;

  void Close() override;

  Result Emit(JEvent& event) override;

  static std::string GetDescription();

  void PrintCollectionTypeTable(void);

  std::vector<std::string_view> getAvailableCategories() const;
  std::size_t getEntries(const std::string& category) const;
  podio::Frame getFrame(const std::string& category, std::size_t index) const;

protected:
  std::unique_ptr<podio::Reader> m_reader;

  std::size_t Nevents_in_file = 0;
  std::size_t Nevents_read    = 0;

  bool m_run_forever       = false;
  bool m_use_event_headers = true;
  bool m_print_type_table  = false;

  // Streaming watch mode (podio:watch_directory): after exhausting the
  // current file, poll a directory for new complete *.edm4hep.root files
  // and continue reading from them instead of finishing. See Emit().
  std::string m_watch_directory;             // empty = off
  std::string m_watch_state_file;            // processed-file ledger, survives rotation
  int m_watch_poll_ms    = 5000;             // min interval between dir scans
  int m_watch_max_files  = 0;                // stop after N extra files; 0 = no limit
  int m_watch_settle_s   = 3;                // ignore files modified more recently than this
  std::set<std::string> m_watch_seen;        // absolute paths already processed
  int m_watch_files_done = 0;
  std::size_t m_watch_frame_offset = 0;      // keeps frame numbers unique across files

  std::chrono::steady_clock::time_point m_watch_last_poll{};

  bool WatchOpenNextFile();                  // returns true when a new file was opened
  void WatchMarkProcessed(const std::string& abs_path); // seen-set + ledger append
  static bool RootFileComplete(const std::string& path); // fEND-vs-size truncation check

  std::shared_ptr<spdlog::logger> m_log;
};

template <> double JEventSourceGeneratorT<JEventSourcePODIO>::CheckOpenable(std::string);
