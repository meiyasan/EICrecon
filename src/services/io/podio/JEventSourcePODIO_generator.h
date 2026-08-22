
#pragma once

#include <JANA/JEventSourceGenerator.h>
#include "JEventSourcePODIO.h"

class JEventSourcePODIO_generator : public JEventSourceGenerator {

  JEventSource* MakeJEventSource(std::string resource_name) override {

    auto* source = new JEventSourcePODIO();
    source->SetTypeName("JEventSourcePODIO");
    source->SetResourceName(resource_name);
    source->SetApplication(mApplication);
    source->SetPluginName(GetPluginName());

    // When the EventBuilder workflow is enabled (-Peventbuilder=true, alias
    // -Peventbld=true), read the input at the Timeslice (time-frame) level so
    // the eventbuilder plugin can build physics events from it; otherwise
    // read it as ordinary physics events.
    const bool use_eventbuilder = mApplication->RegisterParameter<bool>(
        "eventbuilder", false, "Read input as time-frames for the EventBuilder");
    const bool use_eventbld = mApplication->RegisterParameter<bool>(
        "eventbld", false, "Alias of eventbuilder");
    if (use_eventbuilder || use_eventbld) {
      source->SetLevel(JEventLevel::Timeslice);
    } else {
      source->SetLevel(JEventLevel::PhysicsEvent);
    }

    return source;
  }

  double CheckOpenable(std::string resource_name) override {
    // In theory, we should check whether PODIO can open the file and
    // whether it contains an 'events' or 'timeslices' tree. If not, return 0.
    if (resource_name.find(".root") != std::string::npos) {
      return 0.01;
    }
    return 0;
  }

  std::string GetType() const override { return "JEventSourcePodio"; }
};
